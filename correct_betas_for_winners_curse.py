import optparse
import sys
import gzip
import random
import copy
import scipy.optimize
import scipy.stats
import numpy as np

def bail(message):
    sys.stderr.write("%s\n" % (message))
    sys.exit(1)

usage = "usage: correct_betas_for_winners_curse.py --assoc-file <path> --beta-col <int> --se-col <int> --p-threshold <float>"

parser = optparse.OptionParser(usage)
parser.add_option("","--assoc-file",default=[])
parser.add_option("","--beta-col",type="int")
parser.add_option("","--se-col",type="int")
parser.add_option("","--p-col",type="int")
parser.add_option("","--id-col",type="int")
parser.add_option("","--out-delim",default="\t")
parser.add_option("","--p-threshold",type="float")
parser.add_option("","--simulation-reps",type="int",default=100)
parser.add_option("","--simulation-results-out-file")
parser.add_option("","--ci-alpha",type="float",default=0.05)

(options, args) = parser.parse_args()
out_delim = options.out_delim

if len(options.assoc_file) == 0:
    bail(usage)

def open_gz(file):
    if file[-3:] == ".gz":
        return gzip.open(file)
    else:
        return open(file)

def decode(val):
    if type(val) == bytes:
        val = val.decode('utf-8')
    return val

if options.beta_col is None or (options.se_col is None and options.p_col is None) or options.p_threshold is None:
    bail(usage)

class beta_info():
    def __init__(self, beta, se, maf):
        self.beta = beta
        self.se = se

if options.p_threshold <= 0 or options.p_threshold > 1:
    bail("--p-threshold must be between 0 and 1")

threshold_z = scipy.stats.norm.ppf(1 - options.p_threshold / 2)

variant_betas = []

if options.simulation_results_out_file is not None:
    if options.id_col is None:
        bail("When simulating, need to specify --id-col")
    simulation_results_fh = open(options.simulation_results_out_file, 'w')
    simulation_results_fh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("ID", "Beta_MLE", "Beta_MLE_Lower", "Beta_MLE_Upper", "Replicate", "Simulated_Beta", "Simulated_Beta_MSE"))

with open_gz(options.assoc_file) as assoc_fh:
    header = assoc_fh.readline().strip()
    print("%s%sBeta_bias_corrected%sBeta_mle%sBeta_mle_lower%sBeta_mle_upper%sBeta_mse%sBeta_mse_lower%sBeta_mse_upper" % (header, out_delim, out_delim, out_delim, out_delim, out_delim, out_delim, out_delim))
    for line in assoc_fh:
        line = line.strip()
        cols = line.split()
        if options.beta_col > len(cols):
            bail("Not enough columns in %s" % line)
        try:
            beta = float(cols[options.beta_col - 1])
            if options.se_col is not None:
                if options.se_col > len(cols):
                    bail("Not enough columns in %s" % line)
                se = float(cols[options.se_col - 1])
            else:
                if options.p_col > len(cols):
                    bail("Not enough columns in %s" % line)
                p = float(cols[options.p_col - 1])
                z = np.abs(scipy.stats.norm.ppf(p/2))
                if z == 0:
                    continue
                se = abs(beta / z)

        except ValueError:
            continue

        beta_sign = beta / abs(beta)
        beta = abs(beta)

        if beta / se < threshold_z:
            continue

        interval = (0, beta + se)

        beta_true = scipy.optimize.brentq(lambda x: x - beta + se * (scipy.stats.norm.pdf(x / se - threshold_z) - scipy.stats.norm.pdf(-x / se - threshold_z)) / (scipy.stats.norm.cdf(x / se - threshold_z) + scipy.stats.norm.cdf(-x /
 se - threshold_z)), a=interval[0], b=interval[1])
        log_likelihood_func = lambda x: np.log(((1 / se) * scipy.stats.norm.pdf((beta - x)/se))/(scipy.stats.norm.cdf(x/se - threshold_z) + scipy.stats.norm.cdf(-x/se - threshold_z)) * int(np.abs(beta / se) >= threshold_z))

        beta_mle_result = scipy.optimize.minimize(lambda x: -log_likelihood_func(x), beta, bounds=[interval])
        beta_mle = beta_mle_result.x[0]
        beta_mle_fun = -beta_mle_result.fun[0]
        K = se**2 / (se**2 + (beta - beta_mle)**2)
        beta_mse = K * beta + (1 - K) * beta_mle

        #confidence intervals

        if options.ci_alpha <= 0 or options.ci_alpha >= 1:
            bail("--ci-alpha must be between 0 and 1")

        norm_ci_bound_value = abs(scipy.stats.norm.ppf(options.ci_alpha / 2))
        com_ci_lower = beta - norm_ci_bound_value * se
        com_ci_upper = beta + norm_ci_bound_value * se

        mle_ci_bound_value = scipy.stats.chi2.ppf(1 - options.ci_alpha/2, df=1) / 2

        #for x in np.arange(com_ci_lower - 3*se, com_ci_upper + 3 * se, .05):
        #    print("%.3g\t%.3g" % (x, log_likelihood_func(x)))# - mle_ci_bound_value))

        mle_ci_lower = scipy.optimize.brentq(lambda x: log_likelihood_func(x) - (beta_mle_fun - mle_ci_bound_value), a=beta_mle - 10*se, b=beta_mle)
        mle_ci_upper = scipy.optimize.brentq(lambda x: log_likelihood_func(x) - (beta_mle_fun - mle_ci_bound_value), a=beta_mle, b=beta_mle + 10*se)

        K_alpha = se**2 / (se**2 + (com_ci_lower - mle_ci_lower)**2)
        K_1_alpha = se**2 / (se**2 + (com_ci_upper - mle_ci_upper)**2)

        mse_ci_lower = K_alpha * com_ci_lower + (1 - K_alpha) * mle_ci_lower
        mse_ci_upper = K_1_alpha * com_ci_upper + (1 - K_1_alpha) * mle_ci_upper

        if options.simulation_results_out_file is not None:
            chi_statistics = scipy.stats.chi2.rvs(1, size=options.simulation_reps)
            simulation_results = []
            for i in range(len(chi_statistics)):
                se_scale = 4
                while (1):
                    try:
                        beta_lower = scipy.optimize.brentq(lambda x: log_likelihood_func(x) - (beta_mle_fun - chi_statistics[i] / 2), a=beta_mle - 10*se, b=beta_mle)
                        beta_upper = scipy.optimize.brentq(lambda x: log_likelihood_func(x) - (beta_mle_fun - chi_statistics[i] / 2), a=beta_mle, b=beta_mle + 10*se)

                        lower = True if random.random() > 0.5 else False

                        sim_beta = beta_lower if lower else beta_upper

                        chi_quantile = scipy.stats.chi2.sf(chi_statistics[i], df=1)
                        com_matched_beta = beta + (1 if lower else -1) * scipy.stats.norm.ppf(chi_quantile) * se
                        K_sim = se**2 / (se**2 + (com_matched_beta - sim_beta)**2)
                        sim_beta_mse = K_sim * com_matched_beta + (1 - K_sim) * sim_beta
                        sim_beta *= beta_sign
                        sim_beta_mse *= beta_sign
                        simulation_results.append((sim_beta, sim_beta_mse))
                        break
                    except ValueError:
                        se_scale *= 2
                        if se_scale > 100:
                            break

        if beta_sign < 0:
            beta_true = -beta_true
            beta_mle = -beta_mle
            t_mle_ci_lower = -mle_ci_upper
            mle_ci_upper = -mle_ci_lower
            mle_ci_lower = t_mle_ci_lower
            beta_mse = -beta_mse
            t_mse_ci_lower = -mse_ci_upper
            mse_ci_upper = -mse_ci_lower
            mse_ci_lower = t_mse_ci_lower

        print("%s%s%.3g%s%.3g%s%.3g%s%.3g%s%.3g%s%.3g%s%.3g" % (line, out_delim, beta_true, out_delim, beta_mle, out_delim, mle_ci_lower, out_delim, mle_ci_upper, out_delim, beta_mse, out_delim, mse_ci_lower, out_delim, mse_ci_upper
))

        if options.simulation_results_out_file is not None:
            if options.id_col > len(cols):
                bail("Not enough columns in %s" % line)
            row_id = cols[options.id_col -1]
            for i in range(len(simulation_results)):
                simulation_results_fh.write("%s\t%.3g\t%.3g\t%.3g\t%d\t%.3g\t%.3g\n" % (row_id, beta_mle, mle_ci_lower, mle_ci_upper, i+1, simulation_results[i][0], simulation_results[i][1]))


if options.simulation_results_out_file is not None:
    simulation_results_fh.close()
