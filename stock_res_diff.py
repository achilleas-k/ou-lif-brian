from brian import *
import sys
from itertools import combinations, product
import gc

def relerr(act, est):
    return abs(act-est)/act

def find_optimum_noise(archive):
    '''
    Finds the noise amplitude which maximises the difference in output rates
    of two input mu configurations (mu_offs and mu amp).
    '''
    mu_offs_actual = archive['mu_offs_actual']
    mu_offs_est = archive['mu_offs_est']
    mu_amp_actual  = archive['mu_amp_actual']
    mu_amp_est = archive['mu_amp_est']
    sigma_offs_actual = archive['sigma_offs_actual']
    sigma_amp_actual = archive['sigma_amp_actual']
    freq_actual = archive['freq_actual']
    freq_est = archive['freq_est']
    outrates = archive['outrates']
    mo_unique = np.unique(mu_offs_actual)
    ma_unique = np.unique(mu_amp_actual)
    so_unique = np.unique(sigma_offs_actual)
    sa_unique = np.unique(sigma_amp_actual)
    # stick to one frequency
    alldiffs = []
    small_errs = []
    large_errs = []
    for f in np.unique(freq_actual):
        gc.collect()
        f_inds = freq_actual == f
        for ma in ma_unique:
            ma_inds = mu_amp_actual == ma
            mo_pairs = combinations(mo_unique, 2)
            for i, apair in enumerate(mo_pairs):
                small = min(apair)
                large = max(apair)
                # find sigma which produces biggest outrate difference
                small_inds = (mu_offs_actual == small) & ma_inds
                large_inds = (mu_offs_actual == large) & ma_inds
                small_inds = small_inds & f_inds
                large_inds = large_inds & f_inds
                pair_diffs = []
                sigma_confs = product(sa_unique, so_unique)
                for sa, so in sigma_confs:
                    # get the firing rates for each np.unique combo
                    sigma_inds = (sigma_amp_actual == sa) &\
                                 (sigma_offs_actual == so)
                    small_sim_ind = sigma_inds &  small_inds
                    large_sim_ind = sigma_inds &  large_inds
                    small_sim_ind = np.flatnonzero(small_sim_ind)
                    large_sim_ind = np.flatnonzero(large_sim_ind)
                    if len(small_sim_ind) == len(large_sim_ind) > 0:
                        print("Firing rates for ma %f, so %f, sa %f, f %f" % (
                            ma, so, sa, f))
                        print("mo %f: %s" % (small, outrates[small_sim_ind]))
                        print("mo %f: %s" % (large, outrates[large_sim_ind]))
                        rate_diff = max(outrates[large_sim_ind])-\
                                min(outrates[small_sim_ind])
                        pair_diffs.append((sa, so, rate_diff))
                        alldiffs.append(rate_diff)
                        sm_err = relerr(mu_amp_actual[small_sim_ind],
                                        mu_amp_est[small_sim_ind])
                        lg_err = relerr(mu_amp_actual[large_sim_ind],
                                        mu_amp_est[large_sim_ind])
                        small_errs.append(np.count_nonzero(sm_err))
                        large_errs.append(np.count_nonzero(lg_err))
                maxdiff = max(pair_diffs, key=lambda rd: rd[2])
                print("Outrate difference maximized at"
                        " sa %f, so %f. Diff: %f" % (
                    maxdiff[0], maxdiff[1], maxdiff[2]))
                fig = plt.figure()
                fig.clf()
                ax = plt.subplot(111)
                for sa, so, rd in pair_diffs:
                    ax.plot(sa, rd, '^' if rd == maxdiff[2] else 'o',
                            label="%f" % so)
                plt.xlabel('sa')
                plt.ylabel('outrate diff')
                plt.title('mo: [%f, %f], ma: %f, %f' % (small, large, ma, f))
                box=ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.savefig('figures/stoch_res_diff/%f%f%f%f.png' % (
                    small, large, ma, f))
                print("==========================================\n")
        fig = plt.figure()
        fig.clf()
        ax = plt.subplot(111)
        plt.plot(alldiffs, small_errs, 'o', label='small')
        plt.plot(alldiffs, large_errs, 'o', label='large')
        ax.legend(loc='best')
        plt.xlabel('Outrate difference')
        plt.ylabel('Freq. est. error')
        plt.savefig('figures/stoch_res_diff/stochres_vs_fr_est_err_%f.png' % (
            f))
        return


def input_output_ratio(archive):
    '''
    See if there's a correlation between the estimation error and the
    ratio of input to output frequency, for one mu configuration across all
    sigma configurations.
    '''
    mu_peaks_actual = archive['mu_peaks_actual']
    mu_peaks_est = archive['mu_peaks_est']
    mu_offs_actual = archive['mu_offs_actual']
    mu_offs_est = archive['mu_offs_est']
    mu_amp_actual  = archive['mu_amp_actual']
    mu_amp_est = archive['mu_amp_est']
    sigma_offs_actual = archive['sigma_offs_actual']
    sigma_amp_actual = archive['sigma_amp_actual']
    freq_actual = archive['freq_actual']
    freq_est = archive['freq_est']
    outrates = archive['outrates']

    mu_conf_ind = (mu_offs_actual == 0.8) & (mu_amp_actual == 0.8) &\
            (freq_actual == 10)
    print("Using %i simulations" % (count_nonzero(mu_conf_ind)))
    mo = mu_offs_actual[mu_conf_ind]
    ma = mu_amp_actual[mu_conf_ind]
    mp = mu_peaks_actual[mu_conf_ind]
    so = sigma_offs_actual[mu_conf_ind]
    sa = sigma_amp_actual[mu_conf_ind]
    fout = outrates[mu_conf_ind]
    f_err = abs(freq_actual-freq_est)/freq_actual
    mp_err = abs(mu_peaks_actual-mu_peaks_est)/mu_peaks_actual
    ma_err = abs(mu_amp_actual-mu_amp_est)/mu_amp_actual
    mo_err = abs(mu_offs_actual-mu_offs_est)/mu_offs_actual

    plot(mp/fout, mp_err[mu_conf_ind], 'o', label='mp err')
    plot(mo/fout, mo_err[mu_conf_ind], 'o', label='mo err')
    plot(ma/fout, ma_err[mu_conf_ind], 'o', label='ma err')
    legend(loc='best')
    xlabel('outrate')
    ylabel('error')
    show()



if __name__=='__main__':
    filename = sys.argv[1]
    archive = np.load(filename)
    input_output_ratio(archive)
