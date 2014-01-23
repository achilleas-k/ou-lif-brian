import matplotlib
matplotlib.use('Agg') # for plotting remotely or in screen (i.e., no DISPLAY)
from matplotlib import rc
from brian import *

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True) # use latex

plot_save_loc = 'for_nad'

def plot_est_vs_act_err(x_data, y_data, x_label, y_label, figname='',
                                                        plot_zeros=False):
    #rc('text', usetex=True) # use latex
    fontsize = 30
    if not figname:
        figname = '%s_vs_%s' % (x_label, y_label)
    if not len(x_data) or not len(y_data):
        return
    if not plot_zeros:
        nz_i = nonzero(x_data)
        y_data = array(y_data)[nz_i]
        x_data = array(x_data)[nz_i]

    '''
    Compute averages and error bars.
    '''
    x_data = x_data.round(4) # rounding prevents duplicates due to fp precision
    y_data = y_data.round(4)
    x_avgs = []
    y_avgs = []
    y_stds = []
    for ux in unique(x_data):
        x_avgs.append(ux)
        inds = flatnonzero(x_data == ux)
        myavg = mean(y_data[inds])
        mystd = std(y_data[inds])
        y_avgs.append(myavg)
        y_stds.append(mystd)

    x_data = array(x_avgs)
    y_errb = array(y_stds)
    y_data = array(y_avgs)
    oneline = [min(min(y_data), min(x_data)), max(max(y_data), max(x_data))]
    figure(num=1, figsize=(8, 6), dpi=100)
    clf()
    scatter(x_data, y_data, color='black')
    errorbar(x_data, y_data, yerr=y_errb, ecolor='black', fmt=None)
    plot(oneline, oneline, color='gray', linestyle='--', label='1:1')
    xlabel(x_label, fontsize=fontsize)
    xticks(fontsize=fontsize)
    ylabel(y_label, fontsize=fontsize)
    yticks(fontsize=fontsize)
    #legend(loc='best')
    subplots_adjust(left=0.2, bottom=0.15)
    filename_pdf = "%s/%s.pdf" % (plot_save_loc, figname)
    filename_png = "%s/%s.png" % (plot_save_loc, figname)
    savefig(filename_pdf)
    savefig(filename_png)
    print "Figures saved: %s and %s" % (filename_pdf, filename_png)
    clf()

    #err_h, err_b = histogram(relerr, density=True, bins=50)
    #plot(err_b[:-1], err_h)
    #xlabel("RE %s - %s" % (x_label, y_label), fontsize=fontsize)
    #xticks(fontsize=fontsize)
    #axis(xmin=0, xmax=3)
    #filename_pdf = "%s/%s_err.pdf" % (plot_save_loc, figname)
    #filename_png = "%s/%s_err.png" % (plot_save_loc, figname)
    #savefig(filename_pdf)
    #savefig(filename_png)
    #print "Figures saved: %s and %s" % (filename_pdf, filename_png)
    #clf()
    #print "\n\n"
    return



filename = "lotsofdata.npz"
archive = np.load(filename)
mu_offs_est = archive['mu_offs_est']
mu_offs_actual = archive['mu_offs_actual']
mu_amp_est = archive['mu_amp_est']
mu_amp_actual = archive['mu_amp_actual']
mu_peaks_est = archive['mu_peaks_est']
mu_peaks_actual = archive['mu_peaks_actual']
mu_t_est = archive['mu_t_est']
mu_rng_est = archive['mu_rng_est']
sigma_offs_est = archive['sigma_offs_est']
sigma_offs_actual = archive['sigma_offs_actual']
sigma_amp_est = archive['sigma_amp_est']
sigma_amp_actual = archive['sigma_amp_actual']
sigma_peaks_est = archive['sigma_peaks_est']
sigma_peaks_actual = archive['sigma_peaks_actual']
freq_est = archive['freq_est']
freq_actual = archive['freq_actual']
spikecounts = archive['spikecounts']
CV = archive['CV']
mu_peaks_est_bad=archive['mu_peaks_est_bad']

plot_est_vs_act_err(mu_peaks_actual, mu_peaks_est,
        '$\mu_p (mV/ms)$', '$\hat{\mu}_p (mV/ms)$', 'mu_p_est_err',
        plot_zeros=True)
plot_est_vs_act_err(mu_offs_actual, mu_offs_est,
        '$\mu_0 (mV/ms)$', '$\hat{\mu}_0 (mV/ms)$', 'mu_0_est_err',
        plot_zeros=True)
plot_est_vs_act_err(mu_amp_actual, mu_amp_est,
        '$\mu_a (mV/ms)$', '$\hat{\mu}_a (mV/ms)$', 'mu_a_est_err',
        plot_zeros=True)

plot_est_vs_act_err(sigma_peaks_actual, sigma_peaks_est,
        '$\sigma_p (mV/\sqrt{ms})$',
        '$\hat{\sigma}_p (mV/\sqrt{ms})$',
        'sigma_p_est_err', plot_zeros=True)
plot_est_vs_act_err(sigma_offs_actual, sigma_offs_est,
        '$\sigma_0 (mV/\sqrt{ms})$',
        '$\hat{\sigma}_0 (mV/\sqrt{ms})$',
        'sigma_0_est_err', plot_zeros=True)
plot_est_vs_act_err(sigma_amp_actual, sigma_amp_est,
        '$\sigma_a (mV/\sqrt{ms})$',
        '$\hat{\sigma}_a (mV/\sqrt{ms})$',
        'sigma_a_est_err', plot_zeros=True)



