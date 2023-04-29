import optparse
import sys
import gzip
import copy
import scipy.stats
import random
import numpy as np
import statsmodels.api as sm

def bail(message):
    sys.stderr.write("%s\n" % (message))
    sys.exit(1)

usage = "usage: simulate_sample_overlap.py --variant-file <path> --beta-col <int> --se-col <int>"

parser = optparse.OptionParser(usage)
parser.add_option("","--n-case1", type="int")
parser.add_option("","--n-ctrl1",type="int")
parser.add_option("","--n-case12",type="int")
parser.add_option("","--n-ctrl12",type="int")
parser.add_option("","--n-case2", type="int")
parser.add_option("","--n-ctrl2",type="int")
parser.add_option("","--num-reps",type="int", default=100)
parser.add_option("","--num-vars",type="int", default=100)
parser.add_option("","--maf-file")
parser.add_option("","--maf-col",type="int")
parser.add_option("","--top-n1",type="int")
parser.add_option("","--top-n2",type="int")
parser.add_option("","--top-n1-or-top-n2",type="int")
parser.add_option("","--p1",type="float")
parser.add_option("","--p2",type="float")
parser.add_option("","--p1-or-p2",type="float")
(options, args) = parser.parse_args()

def open_gz(file):
    if file[-3:] == ".gz":
        return gzip.open(file)
    else:
        return open(file)

def decode(val):
    if type(val) == bytes:
        val = val.decode('utf-8')
    return val

if options.n_case1 is None or options.n_ctrl1 is None or options.n_case12 is None or options.n_ctrl12 is None or options.n_case2 is None or options.n_ctrl2 is None:
    bail(usage)

class beta_info():
    def __init__(self, beta, se, maf):
        self.beta = beta
        self.se = se
        self.maf = maf

mafs = None
if options.maf_file is not None:
    with open_gz(options.maf_file) as variant_fh:
        header = variant_fh.readline()
        for line in variant_fh:
            cols = line.strip().split()
            if options.maf_col > len(cols):
                bail("Not enough columns in %s" % line)
            try:
                maf = float(cols[options.beta_col - 1])
                mafs.append(maf)
            except ValueError:
                continue

class reg_info():
    def __init__(self, beta1, beta2, p1, p2):
        self.beta1 = beta1
        self.beta2 = beta2
        self.p1 = p1
        self.p2 = p2

results = []
for rep in range(options.num_reps):
    samples = []
    for var in range(options.num_vars):
        if mafs is None:
            maf = scipy.stats.uniform.rvs(size=1)[0]
        else:
            maf = random.choice(mafs)
        case1_genos = scipy.stats.binom.rvs(size=options.n_case1, n=2, p=maf)
        ctrl1_genos = scipy.stats.binom.rvs(size=options.n_ctrl1, n=2, p=maf)
        case12_genos = scipy.stats.binom.rvs(size=options.n_case12, n=2, p=maf)
        ctrl12_genos = scipy.stats.binom.rvs(size=options.n_ctrl12, n=2, p=maf)
        case2_genos = scipy.stats.binom.rvs(size=options.n_case2, n=2, p=maf)
        ctrl2_genos = scipy.stats.binom.rvs(size=options.n_ctrl2, n=2, p=maf)

        study1_y = np.concatenate((np.ones(options.n_case1 + options.n_case12), np.zeros(options.n_ctrl1 + options.n_ctrl12)))
        study2_y = np.concatenate((np.ones(options.n_case2 + options.n_case12), np.zeros(options.n_ctrl2 + options.n_ctrl12)))
        study1_x = np.concatenate((case1_genos, case12_genos, ctrl1_genos, ctrl12_genos))
        study2_x = np.concatenate((case2_genos, case12_genos, ctrl2_genos, ctrl12_genos))

        try:
            model1 = sm.OLS(study1_y, sm.add_constant(study1_x))
            results1 = model1.fit(disp=0)
            beta_tilde1 = results1.params[1]
            t_test_results1 = results1.t_test([[0,1],[1,0]]).summary_frame()
            p1 = np.float64(t_test_results1.loc[['c0'],['P>|%s|' % 't']])[0,0]

            model2 = sm.OLS(study2_y, sm.add_constant(study2_x))
            results2 = model2.fit(disp=0)
            beta_tilde2 = results2.params[1]
            t_test_results2 = results2.t_test([[0,1],[1,0]]).summary_frame()
            p2 = np.float64(t_test_results2.loc[['c0'],['P>|%s|' % 't']])[0,0]

            samples.append(reg_info(beta1=beta_tilde1, beta2=beta_tilde2, p1=p1, p2=p2))
        except IndexError:
            continue

    indices = set(range(len(samples)))

    if options.top_n1 is not None:
        top_indices = sorted(range(len(samples)), key=lambda k: samples[k].p1)[:options.top_n1]
        indices = set([x for x in top_indices if x in indices])
    if options.top_n2 is not None:
        top_indices = sorted(range(len(samples)), key=lambda k: samples[k].p2)[:options.top_n2]
        indices = set([x for x in top_indices if x in indices])
    if options.top_n1_or_top_n2 is not None:
        top_indices1 = sorted(range(len(samples)), key=lambda k: samples[k].p1)[:options.top_n1_or_top_n2]
        top_indices2 = sorted(range(len(samples)), key=lambda k: samples[k].p2)[:options.top_n1_or_top_n2]
        top_indices = set(top_indices1 + top_indices2)
        indices = set([x for x in top_indices if x in indices])
    if options.p1 is not None:
        top_indices = [x for x in range(len(samples)) if samples[x].p1 <= options.p1]
        indices = set([x for x in top_indices if x in indices])
    if options.p2 is not None:
        top_indices = [x for x in range(len(samples)) if samples[x].p2 <= options.p2]
        indices = set([x for x in top_indices if x in indices])
    if options.p1_or_p2 is not None:
        top_indices = [x for x in range(len(samples)) if samples[x].p1 <= options.p1_or_p2 or samples[x].p2 <= options.p1_or_p2]
        indices = set([x for x in top_indices if x in indices])

    num_concordant = 0
    num_study1_larger = 0
    for replicate in indices:
        if samples[replicate].beta1 * samples[replicate].beta2 > 0:
            num_concordant += 1
            if abs(samples[replicate].beta1) > abs(samples[replicate].beta2):
                num_study1_larger += 1

    if len(indices) > 0:
        results.append((float(num_concordant)/len(indices), num_concordant, len(indices), num_study1_larger))

num_concordant = 0
num_concordant2 = 0
frac_concordant = 0
frac_concordant2 = 0
num_study1_larger = 0
num_study1_larger2 = 0
frac_study1_larger = 0
frac_study1_larger2 = 0

for i in range(len(results)):
    cur_results = results[i]
    frac_concordant += cur_results[0]
    frac_concordant2 += cur_results[0] ** 2
    num_concordant += cur_results[1]
    num_concordant2 += cur_results[1] ** 2
    num_study1_larger += cur_results[3]
    num_study1_larger2 += cur_results[3] ** 2
    frac_study1_larger += cur_results[3] / cur_results[1]
    frac_study1_larger2 += (cur_results[3] / cur_results[1]) ** 2

if len(results) > 0:
    print("Concordant: %.2g +/- %.2g (%s/%s)" % (float(frac_concordant)/len(results), np.sqrt(float(frac_concordant2)/len(results) - (float(frac_concordant)/len(results))**2),  float(num_concordant)/len(results), len(indices)))
    print("Study 1 larger and in same direction: %.2g +/- %.2g (%s/%s)" % (float(frac_study1_larger)/len(results), np.sqrt(float(frac_study1_larger2)/len(results) - (float(frac_study1_larger)/len(results))**2), num_study1_larger, nu
m_concordant))
else:
    print("No results passed thresholds")
