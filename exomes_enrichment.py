#Usage: python exomes_enrichment.py
#
#This script applies a one-sided Wilcoxon test to assess whether a gene set has lower than expected p-values
#relative to a matched set of background genes. Selection of background genes can be done at random or by matching to a set of similar genes as defined by metrics in an optional file.
#If effect sizes are specified for both the gene set (expected) and association results (observed), concordance in direction of effect will also be calculated.
#Output is written to stdout
#
#Arguments:
# gene-list-file[path]: the list of genes to test for enrichment. One per row. Must have column for gene ID and optional column for effect size [Required]
# p-value-file[path]: the p-values for the gene set and all other genes in the analysis. Must have columns for gene, p-value, and (optionally) effect size [Required]
# gene-stats-file[path]: an optional file of statistics for each gene to enable better matching to the gene set. Can have as many statistics as desired
# rank-out-file[path]: specify a file to write the ranks of the gene set to
# gene[string]: specify genes to test on the command line (can pass multiple times)
# gene-list-file-gene-col[int]: the 1-based column in the gene-list-file with the gene ID [Default: 1]
# gene-list-file-effect-col[int]: the 1-based column in the gene-list-file with the effect size (beta); range should be (-inf,inf) [Default: None]
# p-value-file-gene-col[int]: the 1-based column in the p-value-file with the gene ID
# p-value-file-p-col[int]: the 1-based column in the p-value-file with the p-value
# p-value-file-effect-col[int]: the 1-based column in the p-value-file with the effect size
# gene-stats-file-gene-col[int]: the 1-based column in the gene-stats-file with the gene ID
# gene-stats-file-stats-col[int]: 1-based columns in the gene-stats-file to use for matching. Specify as many as desired
# --exclude-gene[string]: pass a list of genes to exclude from the gene set
# min-p-for-effect[float]: when effect size concordance is calculated, threshold on p-value for calculation of concordance (e.g. if concordance should only be assessed for p<0.05 associations)
# num-match[int]: the number of background genes to match to each gene in the gene set [Default: 50]
# random-match[bool]: match genes at random [Default: false]
# negative-blacklist-file[path]: a file of gene names to never include in the background set
# warnings-file[path]: write warnings to a file rather than stderr
# debug-level[int]: specify level of output. 0=no log, 1=info, 2=debug, 3=trace
# no-headers[bool]: specify that no files have headers [Default: false]

import numpy
import random
import math
from scipy import spatial
from scipy import stats

import sys

import optparse

def bail(message):
    sys.stderr.write("%s\n" % (message))
    sys.exit(1)

usage = "usage: exomes_enrichment.py --gene-list-file <file> --p-value-file <file> --gene-stats-file <file> [options]"
parser = optparse.OptionParser(usage)

parser.add_option("","--gene-list-file",default=None)
parser.add_option("","--gene-list-file-gene-col",type="int",default=1)
parser.add_option("","--gene-list-file-effect-col",type="int",default=None)
parser.add_option("","--gene",action="append",default=[])
parser.add_option("","--exclude-gene",action="append",default=[])
parser.add_option("","--p-value-file",default=None)
parser.add_option("","--p-value-file-gene-col",type="int",default=None)
parser.add_option("","--p-value-file-p-col",type="int",default=None)
parser.add_option("","--p-value-file-effect-col",type="int",default=None)
parser.add_option("","--min-p-for-effect",type="float",default=1)
parser.add_option("","--num-match",type="int",default=50)
parser.add_option("","--random-match",action="store_true")
parser.add_option("","--negative-blacklist-file",default=None)
parser.add_option("","--gene-stats-file",default=None)
parser.add_option("","--gene-stats-file-gene-col",type="int",default=None)
parser.add_option("","--gene-stats-file-stats-col",type="int",action="append")
parser.add_option("","--warnings-file",default=None)
parser.add_option("","--debug-level",type="int",default=1) #0=no log, 1=info, 2=debug, 3=trace
parser.add_option("","--no-headers",action="store_true")
parser.add_option("","--rank-out-file",default=None)
parser.add_option("","--progressive-mode", action="store_true")
parser.add_option("","--run-norm", action="store_true")
parser.add_option("","--num-lambda", default=100, type="int")
(options, args) = parser.parse_args()

if (options.gene_list_file is None and len(options.gene) == 0) or options.p_value_file is None or options.p_value_file_p_col is None or options.p_value_file_gene_col is None:
    bail(usage)

p_value_file_p_col = options.p_value_file_p_col - 1
p_value_file_gene_col = options.p_value_file_gene_col - 1

gene_list_file_gene_col = options.gene_list_file_gene_col - 1

do_effect_test = False 
p_value_file_effect_col = None
gene_list_file_effect_col = None
if options.gene_list_file_effect_col is not None or options.p_value_file_effect_col is not None:
    if options.gene_list_file_effect_col is None or options.p_value_file_effect_col is None or options.gene_list_file is None:
        bail("If specify effect-col, must do so for both --gene-list-file as well as --p-value-file")
    do_effect_test = True
    p_value_file_effect_col = options.p_value_file_effect_col - 1
    gene_list_file_effect_col = options.gene_list_file_effect_col - 1

if options.gene_stats_file is not None:
    if options.gene_stats_file_gene_col is None or len(options.gene_stats_file_stats_col) == 0:
        bail(usage)
    gene_stats_file_gene_col = options.gene_stats_file_gene_col - 1
    gene_stats_file_stats_cols = list(map(lambda x: x - 1, options.gene_stats_file_stats_col))

warnings_fh = None
if options.warnings_file is not None:
    warnings_fh = open(options.warnings_file, 'w')
else:
    warnings_fh = sys.stderr

def warn(message):
    if warnings_fh is not None:
        warnings_fh.write("%s\n" % message)
        warnings_fh.flush()

NONE=0
INFO=1
DEBUG=2
TRACE=3
def log(message, level=INFO):
    if level <= options.debug_level:
        sys.stderr.write("%s\n" % message)
        sys.stderr.flush()

log("Reading in files...")

p_values_fh = open(options.p_value_file)
if not options.no_headers:
    p_values_header = p_values_fh.readline()
p_values = dict()
effects = dict()
odds_ratio_mode = True
for line in p_values_fh:
    line = line.strip()
    cols = line.split()
    if p_value_file_gene_col < 0 or p_value_file_gene_col >= len(cols):
        bail("--p-value-file-gene-col %s out of bounds" % (p_value_file_gene_col + 1))
    gene = cols[p_value_file_gene_col]
    if p_value_file_p_col < 0 or p_value_file_p_col >= len(cols):
        bail("--p-value-file-p-col %s out of bounds" % (p_value_file_p_col + 1))
    p_value = float(cols[p_value_file_p_col])
    if numpy.isnan(p_value):
        continue
    if gene in p_values:
        bail("Error: gene %s had more than one p-value" % gene)
    p_values[gene] = p_value
    if p_value_file_effect_col is not None and p_value < options.min_p_for_effect:
        if p_value_file_effect_col < 0 or p_value_file_effect_col >= len(cols):
            bail("--p-value-file-effect-col %s out of bounds" % (p_value_file_effect_col + 1))
        effect = float(cols[p_value_file_effect_col])
        effects[gene] = effect
        if effect < 0:
            odds_ratio_mode = False
p_values_fh.close()

positive_genes = {}
positive_gene_effects = dict()
if options.gene_list_file is not None:
    gene_list_fh = open(options.gene_list_file, 'r')
    number = 1
    for line in gene_list_fh:
        line = line.strip()
        cols = line.split()
        if gene_list_file_gene_col < 0 or gene_list_file_gene_col >= len(cols):
            bail("--gene-list-file-gene-col %s out of bounds" % (gene_list_file_gene_col + 1))
        gene = cols[gene_list_file_gene_col]
        if gene not in p_values:
            continue
        positive_genes[gene] = number
        number = number + 1
        if gene_list_file_effect_col is not None:
            if gene_list_file_effect_col < 0 or gene_list_file_effect_col >= len(cols):
                bail("--gene-list-file-effect-col %s out of bounds" % (gene_list_file_effect_col + 1))
            effect = float(cols[gene_list_file_effect_col])
            positive_gene_effects[gene] = effect
            if effect < 0:
                odds_ratio_mode = False
    gene_list_fh.close()

for gene in options.gene:
    if gene in p_values and gene not in options.exclude_gene:
        positive_genes[gene] = number
        number = number + 1

negative_blacklist_genes = set()
if options.negative_blacklist_file is not None:
    blacklist_fh = open(options.negative_blacklist_file, 'r')
    for line in blacklist_fh:
        line = line.strip()
        cols = line.split()
        for gene in cols:
            if gene not in p_values:
                continue
            negative_blacklist_genes.add(gene)
    blacklist_fh.close()

negative_genes = {}
if options.gene_stats_file is not None:

    gene_stats_fh = open(options.gene_stats_file)
    if not options.no_headers:
        gene_stats_header = gene_stats_fh.readline()
    gene_stats_genes = []
    gene_stats_stats = []
    gene_stats = dict()
    for line in gene_stats_fh:
        line = line.strip()
        cols = line.split()
        if gene_stats_file_gene_col < 0 or gene_stats_file_gene_col >= len(cols):
            bail("gene-stats-file-gene-col %s out of bounds" % (gene_stats_file_gene_col+1))
        gene = cols[gene_stats_file_gene_col]
        if gene not in p_values:
            continue
        cur_stats = []
        for gene_stats_file_stats_col in gene_stats_file_stats_cols:
            if gene_stats_file_stats_col < 0 or gene_stats_file_stats_col >= len(cols):
                bail("--gene-stats-file-stats-col %s out of bounds" % (gene_stats_file_stats_col+1))
            cur_stats.append(float(cols[gene_stats_file_stats_col]))
        gene_stats_genes.append(gene)
        gene_stats_stats.append(numpy.array(cur_stats))
        gene_stats[gene] = numpy.array(cur_stats)
    gene_stats_fh.close()

    #normalize
    if len(gene_stats_stats) > 0:
        nums = numpy.zeros(len(gene_stats_stats[0]))
        tots = numpy.zeros(len(gene_stats_stats[0]))
        tots_2 = numpy.zeros(len(gene_stats_stats[0]))
        for i in range(len(gene_stats_stats)):
            tots += gene_stats_stats[i]
            tots_2 += gene_stats_stats[i] ** 2
            nums += numpy.ones(len(gene_stats_stats[i]))
        means = tots / nums
        devs = numpy.sqrt(tots_2 / nums - means ** 2)
        for i in range(len(gene_stats_stats)):
            new_stats = (gene_stats_stats[i] - means) / devs
            gene_stats_stats[i] = new_stats
            gene_stats[gene_stats_genes[i]] = new_stats

    new_positive_genes = {}
    for gene in positive_genes:
        if gene in gene_stats:
            new_positive_genes[gene] = positive_genes[gene]
    positive_genes = new_positive_genes

    kdtree = spatial.cKDTree(numpy.array(gene_stats_stats), leafsize=10)
    for gene in sorted(positive_genes.keys(), key=lambda x: positive_genes[x]):
        results = kdtree.query(gene_stats[gene], options.num_match + 1)
        for i in range(0,min(len(results[1]), len(gene_stats_genes))):
            chosen_gene = gene_stats_genes[results[1][i]]
            if chosen_gene not in positive_genes and chosen_gene not in negative_blacklist_genes:
                if chosen_gene not in negative_genes:
                    negative_genes[chosen_gene] = positive_genes[gene]

elif options.random_match:
    possible_negative_genes = []
    for gene in p_values:
        if gene not in positive_genes and gene not in negative_blacklist_genes:
            possible_negative_genes.append(gene)
    num_to_choose = len(positive_genes) * options.num_match
    if num_to_choose > len(possible_negative_genes):
        num_to_choose = len(possible_negative_genes)

    random.shuffle(possible_negative_genes)
    for i in range(num_to_choose):
        negative_genes[possible_negative_genes[i]] = i / options.num_match

else:
    for gene in p_values:
        if gene not in positive_genes:
            negative_genes[gene] = 0

positive_genes_list = sorted(positive_genes.keys(), key=lambda x: positive_genes[x])
negative_genes_list = sorted(negative_genes.keys(), key=lambda x: negative_genes[x])

#we don't run with only 1 gene so need to seed

cur_positive_genes = set()
if options.progressive_mode:
    slices = range(1, len(positive_genes))
    cur_positive_genes.add(positive_genes_list[0])
    positive_values = [p_values[positive_genes_list[0]]]
    negative_threshold = positive_genes[positive_genes_list[0]]
else:
    slices = [len(positive_genes)]
    for gene in positive_genes_list:
        cur_positive_genes.add(gene)
    positive_values = [p_values[x] for x in positive_genes_list]
    negative_threshold = None

positive_genes_info = dict()

cur_negative_genes = set()
negative_values = []

negative_genes_index = 0

positive_ranks = []
R1 = 0
n1 = 0
n2 = 0

for i in slices:

    if options.progressive_mode:
        cur_positive_genes.add(positive_genes_list[i])
        positive_values.append(p_values[positive_genes_list[i]])
        negative_threshold = positive_genes[positive_genes_list[i]]
    while negative_genes_index < len(negative_genes_list) and (negative_threshold is None or negative_genes[negative_genes_list[negative_genes_index]] <= negative_threshold):
        cur_negative_genes.add(negative_genes_list[negative_genes_index])
        negative_values.append(p_values[negative_genes_list[negative_genes_index]])
        negative_genes_index += 1

    result = stats.mannwhitneyu(positive_values, negative_values)

    rank_out_fh = None
    if options.rank_out_file and not options.progressive_mode:
        rank_out_fh = open(options.rank_out_file, 'w')
        rank_out_fh.write("Gene\tP-value\tRank\tPercentile\tPositive\n")

    rank = 1
    num_neg = 0
    num_pos = 0
    avg_log_p_value = 0
    for gene in sorted(p_values.keys(), key=lambda x: p_values[x]):
        percentile = 1 - (float(rank) / (len(cur_positive_genes) + len(cur_negative_genes)))
        if gene in cur_positive_genes:
            if rank_out_fh is not None:
                rank_out_fh.write("%s\t%.3g\t%d\t%.3g\t%d\n" % (gene, p_values[gene], rank, percentile, 1))
            positive_genes_info[gene] = (p_values[gene], rank)
            R1 += rank
            n1 += 1
            rank += 1
            num_pos += 1
        elif gene in cur_negative_genes:
            if rank_out_fh is not None:
                rank_out_fh.write("%s\t%.3g\t%d\t%.3g\t%d\n" % (gene, p_values[gene], rank, percentile, 0))
            n2 += 1
            rank += 1
            avg_log_p_value += numpy.log(p_values[gene])
            num_neg += 1
    avg_log_p_value /= num_neg
    fold_frac = 0.1

    positive_p_values = [x[1][0] for x in positive_genes_info.items()]

    fold_threshold = sorted([p_values[x] for x in p_values])[int(fold_frac * len(p_values))]
    negative_fold = len([p_values[x] for x in p_values if x in cur_negative_genes and p_values[x] < fold_threshold]) / float(num_neg)
    positive_fold = len([x for x in positive_p_values if x < fold_threshold]) / float(num_pos)
    
    if rank_out_fh is not None:
        rank_out_fh.close()

    U1 = R1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    sigmau = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)

    Z = (U1 - mu) / sigmau
    from scipy import stats
    norm_pvalue = stats.norm.cdf(Z)

    if do_effect_test:
        num_same_direction = 0
        num_diff_direction = 0
        for gene in positive_gene_effects:
            if gene in effects:
                ref_effect = positive_gene_effects[gene]
                obs_effect = effects[gene]
                if odds_ratio_mode:
                    ref_effect = ref_effect - 1
                    obs_effect = obs_effect - 1
                if ref_effect * obs_effect > 0:
                    num_same_direction += 1
                elif ref_effect * obs_effect < 0:
                    num_diff_direction += 1
        binom_pvalue = stats.binom_test(num_same_direction, n=num_same_direction + num_diff_direction, p=0.5) / 2
        if num_same_direction < num_diff_direction:
            binom_pvalue = 1 - binom_pvalue

    result_pvalue = result.pvalue
    if abs(1 - norm_pvalue - result.pvalue) < abs(norm_pvalue - result.pvalue):
        result_pvalue = 1 - result.pvalue

    # W(lambda)=number{p : p > lambda}
    # pi(lambda) = W(lambda) / (1 - lambda) m
    tot_W = 0
    tot_l = 0
    j_star = None
    intervals = numpy.linspace(0,1,options.num_lambda)
    for j in range(len(intervals))[:-1]:
        l = intervals[j]
        l_p_1 = intervals[j+1]
        b_j = float(sum([x >= l for x in positive_p_values]))
        c_j = float(sum([x >= l and x < l_p_1 for x in positive_p_values]))
        if j_star is None and c_j <= b_j / (len(intervals) - j):
            j_star = j
        tot_W += b_j / ((1 - l) * len(positive_p_values))
        
    prop_true = 1 - (tot_W / (len(intervals) - j_star))

    mean_log_p_ratio = numpy.log(numpy.mean([numpy.log(x) for x in positive_p_values if x > 0]) / avg_log_p_value)

    fold_enrichment = numpy.inf
    if negative_fold > 0:
        fold_enrichment = positive_fold / negative_fold
    if options.progressive_mode:
        if i == slices[0]:
            print("Rank\tP\tProp_true\tLog_ratio\tFold")
        print("%d\t%.3g\t%.3g\t%.3g\t%.3g" % (i+1, result_pvalue, prop_true, mean_log_p_ratio, fold_enrichment))
    else:
        print("Number positive: %d" % len(positive_values))
        print("Number negative: %d" % len(negative_values))
        if options.gene_stats_file:
            print("Matched on variables in %s" % options.gene_stats_file)
        print("Significance: p=%.4g" % result_pvalue)
        print("Norm significance: p=%.4g" % norm_pvalue)
        print("Prop true: %.4g" % prop_true)
        print("Log mean log(p) ratio: %.4g" % mean_log_p_ratio)
        print("Fold enrichment (top %s): : %.4g" % (fold_frac, fold_enrichment))
        if do_effect_test:
            print("Directional sign test: p=%.4g (%d vs. %d)" % (binom_pvalue, num_same_direction, num_diff_direction))
        header = "%s\t%s\t%s\t%s\t%s" % ("Gene", "P-value", "Null_P-value", "Rank", "Percentile")
        if do_effect_test:
            header = "%s\t%s\t%s" % (header, "Reference_Sign", "Observed_Sign")
        print(header)
        pos = 1
        for x in sorted(positive_genes_info, key=lambda x: positive_genes_info[x][0]):
            output = "%s\t%.3g\t%.3g\t%d\t%.3g" % (x, positive_genes_info[x][0], float(pos) / len(positive_values), positive_genes_info[x][1], 1-float(positive_genes_info[x][1]) / (len(positive_values) + len(negative_values)))
            if do_effect_test:
                output = "%s\t%d\t%d" % (output, positive_gene_effects[x] / abs(positive_gene_effects[x]) if x in positive_gene_effects and positive_gene_effects[x] != 0 else 0, effects[x] / abs(effects[x]) if x in effects and effec
ts[x] != 0 else 0)
            print(output)
            pos += 1
