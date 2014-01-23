'''
Plot all the figures for the 2013 Brain Research paper as they will appear
in the manuscript.
'''
from __future__ import division
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc
from brian import *
import sys, os


#rc('font',**{'family':'sans-serif','sans-serif':['Bitstream Vera Sans']})
#rc('text', usetex=True) # use latex
plot_save_loc = 'brain_research_figures'

def plot_methodology():
    '''
    Generates Figure 1 in the paper which outlines the estimation methodology
    '''
    defaultclock.dt = dt = 0.1*ms
    freq = 10*Hz
    period = 1/freq
    duration = 4*period
    mu_offs = 3*mV/ms
    mu_amp = 1*mV/ms
    sigma_offs = 1*mV/sqrt(ms)
    sigma_amp = 1*mV/sqrt(ms)
    T = frange(0*second, duration, 0.1*ms)
    signal = mu_amp*sin(freq*T*2*pi)+mu_offs
    ncoef = sigma_amp*sin(freq*T*2*pi)+sigma_offs
    noise = normal(0, 1, len(T))*ncoef
    # let's get random samples from the original sine wave
    # the larger the value, the higher the probability of getting a sample
    probabilities = (signal-min(signal))/max(signal)
    probabilities *= 0.3   # so that peaks are not p=1
    T_est = []
    # the following isn't the simplest way to do it, but allows us to check
    # for refractory period
    for ti, p in zip(T, probabilities):
        if rand() < p:
            if not T_est or ti-T_est[-1] >= 0.002:
                T_est.append(ti)
    T_est = array(T_est)
    est = mu_amp*sin(freq*T_est*2*pi)+mu_offs
    rel_errors = random(len(est))*0.2-0.1
    est += est*rel_errors
    duraticks = frange(period, duration, period)
    duraticklabels = ['$%i f^{-1}$' % (i+1) if i else '$ f^{-1} $'
                                            for i in range(len(duraticks))]
    periodticks = [0, float(period)/2, float(period)]
    periodticklabels = ['0', '$0.5 f^{-1}$', '$f^{-1}$']
    muticks = [float(mu_offs-mu_amp), float(mu_offs), float(mu_offs+mu_amp)]
    muticklabels = ['$\mu_0-\mu_a$', '$\mu_0$', '$\mu_p$']
    periodspikes = T_est % period
    binedges = linspace(0*second, period, 11)
    avgest = []
    for l, r in zip(binedges[:-1], binedges[1:]):
        binspike_inds = bitwise_and(periodspikes > l, periodspikes < r)
        binmu = est[binspike_inds]
        avgest.append(mean(binmu) if len(binmu) else 0)
    avgest = array(avgest)
    axeslimits_full = [0, float(duration),
            (mu_offs-mu_amp)*0.95, (mu_offs+mu_amp)*1.05]
    axeslimits_period = [0, float(period),
            (mu_offs-mu_amp)*0.95, (mu_offs+mu_amp)*1.05]
    titalign = ' '*110 # spaces to align title to the left
    # start plotting
    # build a rectangle in axes coords
    figure(1, figsize=(8, 6), dpi=100)
    clf()
    subplot(5, 1, 1)
    plot(T, signal+noise, 'k-')
    plot(T, signal, color='grey', linestyle='--')
    title(r'A%s' % (titalign))
    #title('Original signal (with noise)')
    xticks(duraticks, duraticklabels)
    yticks(muticks, muticklabels)
    axis(axeslimits_full)
    subplot(5, 1, 2)
    stem(T_est, ones(len(T_est))*(mu_offs+mu_amp/2), linefmt='k-', basefmt='k-',
            markerfmt='w-', bottom=mu_offs)
    plot(T, signal, color='grey', linestyle='--')
    title(r'B%s' % titalign)
    #title('Fired spikes')
    xticks(duraticks, duraticklabels)
    yticks(muticks, muticklabels)
    #yticks([])
    axis(axeslimits_full)
    subplot(5, 1, 3)
    plot(T_est, est, 'k.')
    plot(T, signal, color='grey', linestyle='--')
    title(r'C%s' % titalign)
    #title('Estimated $\mu$ values $\hat{\mu}_i$')
    xticks(duraticks, duraticklabels)
    yticks(muticks, muticklabels)
    axis(axeslimits_full)
    subplot(5,1,4)
    plot(periodspikes, est, 'k.')
    plot(T, signal, color='grey', linestyle='--')
    title(r'D%s' % titalign)
    #title('Estimated $\mu$ values, aligned to one period $\hat{\mu}_i$')
    xticks(periodticks, periodticklabels)
    yticks(muticks, muticklabels)
    axis(axeslimits_period)
    subplot(5,1,5)
    plot(binedges[:-1]+0.5*period/10, avgest, color='black', linestyle='-',
            marker='.', markersize=10)
    plot(T, signal, color='grey', linestyle='--')
    title(r'E%s' % titalign)
    #title('Binned average of estimated $\mu$ values $\hat{\mu}_b$')
    xticks(periodticks, periodticklabels)
    yticks(muticks, muticklabels)
    axis(axeslimits_period)
    subplots_adjust(hspace=0.9)
    savefig('%s/method.svg' % (plot_save_loc))
    savefig('%s/method.png' % (plot_save_loc))
    savefig('%s/method.pdf' % (plot_save_loc))
    savefig('%s/method.eps' % (plot_save_loc))
    print("Saved method figure.")
    return

def plot_est_vs_act_err(x_data, y_data, x_label, y_label, figname='',
                                axis_limits=None, plot_zeros=False,
                                size=(800,600), plot_title='',
                                ticks=None):
    '''
    Used to draw the seven parameter estimation plots:
    - Frequency
    - mu_p, mu_o, mu_a
    - sigma_p, sigma_o, sigma_a
    '''
    fontsize = 20
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
    dpi = 100
    size = (size[0]/dpi, size[1]/dpi)
    figure(figsize=size, dpi=dpi)
    scatter(x_data, y_data, color='black')
    errorbar(x_data, y_data, yerr=y_errb, ecolor='black', fmt=None)
    plot(oneline, oneline, color='gray', linestyle='--', label='1:1')
    xlabel(x_label, fontsize=fontsize)
    if ticks is None:
        xticks(fontsize=fontsize)
        yticks(fontsize=fontsize)
    else:
        xticks(ticks[0], ticks[1], fontsize=fontsize)
        yticks(ticks[0], ticks[1], fontsize=fontsize)
    ylabel(y_label, fontsize=fontsize)
    title(plot_title, fontsize=fontsize)
    subplots_adjust(left=0.175, bottom=0.25, right=0.85)
    if axis_limits is not None:
        axis(axis_limits)
    savetypes = ["pdf", "png", "eps", "svg"]
    for ft in savetypes:
        filename = "%s/%s.%s" % (plot_save_loc, figname, ft)
        savefig(filename)
        print("Figures saved: %s" % (filename))
    clf()
    return

def mu_offs_amp_fee(mu_offs_actual, mu_amp_actual, freq_actual,
        freq_est):
    '''
    mu offs vs mu amp for non-zero frequency estimation errors
    '''
    fe_err = abs(freq_actual-freq_est)/freq_actual
    figure(num=1, figsize=(8, 6), dpi=100)
    clf()
    nzerr = flatnonzero(fe_err)
    scatter(mu_offs_actual, mu_amp_actual, c='black', marker='.', s=10,
            label=r'$ \varepsilon_f  = 0 $')
    scatter(array(mu_offs_actual)[nzerr], array(mu_amp_actual)[nzerr],
            c='black', marker='^', s=100, label=r'$ \varepsilon_f > 0 $')
    xlabel('$\mu_0 \mathrm{(mV/ms)}$', size=15)
    ylabel('$\mu_a \mathrm{(mV/ms)}$', size=15)
    xticks(fontsize=15)
    yticks(fontsize=15)
    legend()
    title("Non-zero frequency estimation errors", size=15)
    savefig('%s/mu_for_nz_freq_err.svg' % plot_save_loc)
    savefig('%s/mu_for_nz_freq_err.eps' % plot_save_loc)
    savefig('%s/mu_for_nz_freq_err.png' % plot_save_loc)
    savefig('%s/mu_for_nz_freq_err.pdf' % plot_save_loc)

def plot_cv_errf(CV, freq_est, freq_actual):
    '''
    Plots frequency estimation error vs CV and conditional probability of
    error on CV.
    '''
    ferr = abs(freq_est-freq_actual)/freq_actual
    figure(num=1, figsize=(8, 6), dpi=100)
    clf()
    fontsize = 15
    subplot(2,1,1)
    scatter(CV, ferr, c='k')
    ylabel(r'$ \varepsilon_f $', size=fontsize)
    xticks(size=fontsize)
    yticks(frange(2, 10, 2), size=fontsize)
    axis([0, 2.5, 0, 10])
    subplot(2,1,2)
    cv_bins = frange(0, 2.5, 0.05)
    CVh, ignore = histogram(CV, bins=cv_bins)
    CVe, ignore = histogram(CV[ferr>0], bins=cv_bins)
    condProb = ([1.0*cve/cvh if cvh else 0 for cve, cvh in zip(CVe, CVh)])
    plot(cv_bins[:-1], condProb, 'k')
    xlabel('CV', size=fontsize)
    xticks(size=fontsize)
    yticks(frange(0.2, 0.6, 0.2), size=fontsize)
    savefig('%s/CV_freq_est_err.eps' % (plot_save_loc))
    savefig('%s/CV_freq_est_err.svg' % (plot_save_loc))
    savefig('%s/CV_freq_est_err.png' % (plot_save_loc))
    savefig('%s/CV_freq_est_err.pdf' % (plot_save_loc))
    print("Saved CV freq est err")
    return

def plot_example_results(mu_t_est, mu_offs, mu_amp,
                        sigma_offs, sigma_amp, freq, spiketimes):
    '''
    Plots three figures that display the methodology on example data.
    Also plot a fourth, three-panel figure as an alternative.
    '''
    # find simulation index
    index = (mu_offs == 1.6) & (mu_amp == 0.4) &\
            (sigma_offs == 0.4) & (sigma_amp == 0.4) & (freq == 5)
    if len(index) > 1:
        print('Multiple simulations found for example data.')
        index = index[0]
    mo = mu_offs[index]
    ma = mu_amp[index]
    so = sigma_offs[index]
    sa = sigma_amp[index]
    f  = freq[index]
    spikes = spiketimes[index]
    mt = mu_t_est[index]
    duration = 5*second
    t = arange(0*second, duration, 0.1*ms)
    theoretical = sin(f*2*pi*t)*ma+mo
    # individual estimates
    figure(1)
    plot(spikes, mt, '-o', color='black', label=r'$\hat\mu$')
    plot(t, theoretical, '--', color='grey', label=r'$\mu(t)$')
    legend(loc='best')
    xlabel(r't (sec)')
    ylabel(r'$\mu (mV/st)$')
    axis(xmax=1)
    savefig('%s/mu_t_est_actual.eps' % (plot_save_loc))
    savefig('%s/mu_t_est_actual.png' % (plot_save_loc))
    savefig('%s/mu_t_est_actual.svg' % (plot_save_loc))
    savefig('%s/mu_t_est_actual.pdf' % (plot_save_loc))
    print('Saved mu_t_est_actual.')
    # grouped estimates
    figure(2)
    period = 1/f
    periodspikes = spikes % period
    plot(periodspikes, mt, 'o', color='black', label=r'$\hat\mu$')
    plot(t, theoretical, '--', color='grey', label=r'$\mu(t)$')
    legend(loc='best')
    xlabel(r't (sec)')
    ylabel(r'$\mu (mV/st)$')
    axis(xmin=0, xmax=period)
    savefig('%s/mu_t_est_actual_period.eps' % (plot_save_loc))
    savefig('%s/mu_t_est_actual_period.png' % (plot_save_loc))
    savefig('%s/mu_t_est_actual_period.svg' % (plot_save_loc))
    savefig('%s/mu_t_est_actual_period.pdf' % (plot_save_loc))
    print('Saved mu_t_est_actual_period.')
    # binned averages
    figure(3)
    binedges = linspace(0*second, period, 11)
    avgest = []
    for l, r in zip(binedges[:-1], binedges[1:]):
        binspike_inds = (periodspikes > l) & (periodspikes < r)
        binmu = est[binspike_inds]
        avgest.append(mean(binmu) if len(binmu) else 0)
    avgest = array(avgest)
    plot(binedges[:-1]+0.5*period/10, avgest, color='black', linestyle='-',
            marker='.', markersize=10)
    plot(t, theoretical, '--', color='grey', label=r'$\mu(t)$')
    savefig('%s/mu_t_est_actual_binned.eps' % (plot_save_loc))
    savefig('%s/mu_t_est_actual_binned.png' % (plot_save_loc))
    savefig('%s/mu_t_est_actual_binned.pdf' % (plot_save_loc))
    savefig('%s/mu_t_est_actual_binned.svg' % (plot_save_loc))
    print('Saved mu_t_est_actual_binned.')
    return

def discrete_mu_recon(spikes, V=None, dt=0.1*ms):
    '''
    Estimates the value of mu for each interspike interval based on the
    mu estimator from Bibbona et al., 2008.
    '''
    recon = []
    st = append(0, spikes)
    for tp, t in zip(st[:-1], st[1:]):
        t *= second
        tp *= second
        isi = t-tp
        K = round(isi/dt)
        Vi = V[(tp/dt+1):(t/dt)]
        me = V_th / (K*tau*(1-exp(-dt/tau))) + 1/(K*tau) * sum(Vi)*volt
        recon.append(me)
    return array(recon)

def discrete_mu_recon_avg(spikes, V, freq, nbins, dt=0.1*ms):
    '''
    Estimates the value of mu for each interspike interval based on the mu
    estimator from Bibbona et al., 2008, folds the estimated values
    based on period using the frequency `freq` and averages the estimated
    values in `nbins` bins.

    TODO: Try a moving average instead of binning.
    '''
    recon = discrete_mu_recon(spikes, V, dt)
    if is_dimensionless(freq):
        freq *= Hz
    period = 1/freq
    periodspikes = spikes % period
    binedges = linspace(0*second, period, nbins+1)
    avgmu = []
    for l, r in zip(binedges[:-1], binedges[1:]):
        binspike_inds = bitwise_and(periodspikes > l, periodspikes < r)
        binmu = recon[binspike_inds]
        avgmu.append(mean(binmu) if len(binmu) else 0)
    return array(avgmu)



if __name__=='__main__':
    # Load the data
    filename = 'lotsofdata.npz'
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
    mu_peaks_est_bad = archive['mu_peaks_est_bad']
    outrates = archive['outrates']
    spiketimes = archive['spiketimes']

    # Start plotting
    # Figure 1
    #plot_methodology()
    # Figure 2
    #plot_est_vs_act_err(freq_actual, freq_est, '$f \mathrm{(Hz)}$',
    #        '$\hat{f} \mathrm{(Hz)}$',
    #    'frequency_estimation', plot_zeros=True, axis_limits=[0, 25, 0, 25])
    # Figure 3
    #mu_offs_amp_fee(mu_offs_actual, mu_amp_actual, freq_actual, freq_est)
    # Figure 4
    plot_cv_errf(CV, freq_est, freq_actual)
    # Figure 5, 6, 7
    #plot_example_results(mu_t_est, mu_offs_actual, mu_amp_actual,
    #            sigma_offs_actual, sigma_amp_actual, freq_actual, spiketimes)
    # Figure 8A
    #plot_est_vs_act_err(mu_peaks_actual, mu_peaks_est,
    #        '$\mu_p \mathrm{(mV/ms)}$',
    #        '$\hat{\mu}_p \mathrm{(mV/ms)}$', 'mu_peaks_estimation',
    #        plot_zeros=True,
    #        axis_limits=[0, 4.5, 0, 4.5], size=(1200, 450),
    #        plot_title='A'+' '*110,
    #        ticks=(frange(0, 4.5, 0.5), [0, '', 1, '', 2, '', 3, '', 4, '']),
    #        )
    ## Figure 8B
    #plot_est_vs_act_err(mu_offs_actual, mu_offs_est,
    #        '$\mu_0 \mathrm{(mV/ms)}$',
    #        '$\hat{\mu}_0 \mathrm{(mV/ms)}$', 'mu_offs_estimation',
    #        plot_zeros=True,
    #        axis_limits=[0, 2.5, 0, 2.5],
    #        size=(600, 450),
    #        plot_title='B'+' '*65)
    ## Figure 8C
    #plot_est_vs_act_err(mu_amp_actual, mu_amp_est, '$\mu_a \mathrm{(mV/ms)}$',
    #        '$\hat{\mu}_a \mathrm{(mV/ms)}$',
    #        'mu_amp_estimation', plot_zeros=True,
    #        axis_limits=[0, 2.5, 0, 2.5],
    #        size=(600, 450),
    #        plot_title='C'+' '*65)



    ## Figure 9A
    #plot_est_vs_act_err(sigma_peaks_actual, sigma_peaks_est,
    #        '$\sigma_p \mathrm{(mV/\sqrt{ms})}$',
    #        '$\hat{\sigma}_p \mathrm{(mV/\sqrt{ms})}$',
    #        'sigma_peaks_estimation',
    #        plot_zeros=True, axis_limits=[-0.2, 2.5, -0.2, 2.5],
    #        size=(1200, 450), plot_title='A'+' '*110,
    #        ticks=(frange(0, 2.5, 0.25), [0, '', 0.5, '', 1, '', 1.5, '', 2]),
    #        )

    ## Figure 9B
    #plot_est_vs_act_err(sigma_offs_actual, sigma_offs_est,
    #        '$\sigma_0 \mathrm{(mV/\sqrt{ms})}$',
    #        '$\hat{\sigma}_0 \mathrm{(mV/\sqrt{ms})}$',
    #        'sigma_offs_estimation',
    #        plot_zeros=True, axis_limits=[-0.2, 1.2, -0.2, 1.2],
    #        size=(600, 450),
    #        plot_title='B'+' '*65,
    #        ticks=(frange(0, 1, 0.25), [0, '', 0.5, '', 1]),
    #        )
    ## Figure 9C
    #plot_est_vs_act_err(sigma_amp_actual, sigma_amp_est,
    #        '$\sigma_a \mathrm{(mV/\sqrt{ms})}$',
    #        '$\hat{\sigma}_a \mathrm{(mV/\sqrt{ms})}$',
    #        'sigma_amp_estimation',
    #        plot_zeros=True, axis_limits=[-0.2, 1.2, -0.2, 1.2],
    #        size=(600, 450),
    #        plot_title='C'+' '*65,
    #        ticks=(frange(0, 1, 0.25), [0, '', 0.5, '', 1]),
    #        )





