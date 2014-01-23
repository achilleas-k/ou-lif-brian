import matplotlib
matplotlib.use('Agg') # for plotting remotely or in screen (i.e., no DISPLAY)
from matplotlib import rc
from brian import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import sys
import gc
import itertools
from time import time
import neurotools as nt


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True) # use latex
'''
Intrinsic neuron parameters which are assumed to be known
'''
V_th = 10*mV
tau = 10*ms
t_refr = 0*ms
v_reset = 0*mV
V0 = 0*mV
dt = 0.1*ms

plot_save_loc = "brain_research_figures"

def plot_est_vs_act(x_data, y_data, x_label, y_label, figname='',
                                                        plot_zeros=False):
    fontsize = 15
    if not figname:
        figname = '%s_vs_%s' % (x_label, y_label)
    if not len(x_data) or not len(y_data):
        return
    if not plot_zeros:
        nz_i = nonzero(x_data)
        y_data = array(y_data)[nz_i]
        x_data = array(x_data)[nz_i]
    good_thr = 0.1
    x_data = array(x_data)
    y_data = array(y_data)
    relerr = abs(x_data-y_data)/x_data
    good_ratio = 1.0*count_nonzero(relerr < good_thr)/len(relerr)
    oneline = [min(min(y_data), min(x_data)), max(max(y_data), max(x_data))]
    figure(num=1, figsize=(8, 6), dpi=100)
    clf()
    scatter(x_data, y_data, c='black', label=x_label)
    plot(oneline, oneline, color='gray', linestyle='--', label='1:1')
    xlabel(x_label, fontsize=fontsize)
    xticks(fontsize=fontsize)
    ylabel(y_label, fontsize=fontsize)
    yticks(fontsize=fontsize)
    #legend(loc='best')
    print('%s vs %s.\n\t\t%f%% have an error < %0.2f '
            '\n\t\tMean error: %f' % (
        y_label, x_label,
        good_ratio*100,
        good_thr,
        mean(relerr)))
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

def plot_est_vs_act_err(x_data, y_data, x_label, y_label, figname='',
                                                        plot_zeros=False):
    #rc('text', usetex=True) # use latex
    fontsize = 15
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
    x_data = x_data.round(2) # rounding prevents duplicates due to fp precision
    y_data = y_data.round(2)
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

def psd_est(st, duration):
    '''
    Estimate the frequency based on the fourier transform of the
    autocorrelogram (power spectrum) of the spike train
    '''
    width = min(duration, 2*second)
    bin = 10*ms # best so far = 10*ms
    acorr = autocorrelogram(st, width=width, bin=bin, T=duration)
    acorr = acorr[len(acorr)/2:] # acorr is symmetric
    psd = fft(acorr)
    cutfreq = len(psd)/2
    psd[cutfreq:] = 0 # low pass filter
    psd[0] = 0
    freq_est_i = argmax(psd)
    freq_est = (1.0*freq_est_i/len(psd))/bin
    return freq_est

def find_burst_length(spiketrain):
    '''
    Find the burst length, i.e., the mean consecutive number of spikes with
    an interval < mean ISI.
    '''
    ISIs = diff(spiketrain)
    mISI = mean(ISIs)
    sub_mISI = flatnonzero(ISIs < mISI) # new flatnonzero
    cons = 0
    ncons = []
    for ii in diff(sub_mISI):
        if ii == 1:
            cons += 1
        else:
            ncons.append(cons)
            cons = 0
    if len(ncons) < 2:
        return 0
    return int(mean(ncons)+1)

def scaled_period_spikes(spiketrain, freq):
    if is_dimensionless(freq):
        freq = freq*Hz
    period = 1.0/freq
    period_spikes = spiketrain % period
    nperiods = int(round(spiketrain[-1]/period))
    period_binned = nt.times_to_bin(period_spikes)/nperiods
    t = arange(0*second, period, 1*ms)
    period_binned = append(period_binned, zeros(len(t)-len(period_binned)))
    scaled_spikes = period_binned/(sin(t*2*pi*freq)+1.1)
    return scaled_spikes

def discrete_mu_recon(spikes, V=None, dt=0.1*ms):
    '''
    Estimates the value of mu for each interspike interval based on the
    mu estimator from Bibbona et al., 2008.
    '''
    isi = diff(spikes)
    isi = insert(isi, 0, spikes[0])
    recon = []
    for T, sp in zip(isi, spikes):
        K = T/dt
        if V is not None:
            iii = (arange((sp-T)*second, (sp)*second, dt)/dt).astype(int)
            memsum = sum(V[iii])*volt
        else:
            memsum = sum(arange(0*second, T*second, dt))*volt
        m_est = V_th/(tau*K*(1-exp(-dt/tau))) + memsum/(tau*K)
        recon.append(m_est)


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
    bins = linspace(0*second, period, nbins)
    avgmu = []
    for l, r in zip(bins[:-1], bins[1:]):
        binspike_inds = bitwise_and(periodspikes > l, periodspikes < r)
        binmu = recon[binspike_inds]
        avgmu.append(mean(binmu) if len(binmu) else 0)
    return array(avgmu)

def slope_mu_recon(spikes, V, freq, binwidth, w, dt=0.1*ms):
    msl, slopes = nt.firing_slope(V, spikes, w, dt)
    if not len(slopes):
        return zeros(1/(freq*binwidth))
    if is_dimensionless(binwidth):
        binwidth *= second
    if is_dimensionless(freq):
        freq *= Hz
    period = 1/freq
    periodspikes = spikes % period
    bins = frange(0*second, period, binwidth)
    avgslopes = []
    for l, r in zip(bins[:-1], bins[1:]):
        binspike_inds = bitwise_and(periodspikes > l, periodspikes < r)
        binslopes = slopes[binspike_inds]
        avgslopes.append(mean(binslopes) if len(binslopes) else 0)
    return array(avgslopes)

def is_super(spikes, duration, freq, tau): # version 1
    '''
    Estimates whether the mu offs value that caused `spikes` is super- or sub-
    threshold.
    '''
    period = 1/freq
    periodspikes = spikes % period
    highend = 0.30*period
    lowstart = 0.65*period
    nhighspikes = count_nonzero(periodspikes < highend)
    midspikes = count_nonzero(bitwise_and(highend < periodspikes,
                periodspikes < lowstart))
    lowspikes = count_nonzero(periodspikes > lowstart)
    nperiods = duration/period
    if lowspikes:
        return 1
    elif midspikes:
        firstspike = sorted(periodspikes)[0]*second
        if firstspike > tau: # threshold should be based on vth/tau
            return 0
        else:
            return 1
    else:
        return 0

if __name__ == '__main__':
    warn('WARNING: Unmaintained script. Use InpEstimationMP instead')
    sys.exit(1)
    data_dir = sys.argv[1]
    if data_dir.endswith('.data'):
        data_dir = data_dir.replace('.data', '')
    elif data_dir.endswith('.data/'):
        data_dir = data_dir.replace('.data/', '')
    data = DataManager(data_dir)
    print("\n")
    print "Working with %s.data directory\n" %\
                data_dir

    numsims = data.itemcount()
    print "Total number of simulations: %i" % numsims

    '''
    mu lists
    '''
    mu_offs_est = []
    mu_offs_actual = []

    mu_amp_est = []
    mu_amp_actual = []

    mu_peaks_est = []
    mu_peaks_actual = []

    mu_t_est = []
    mu_rng_est = []

    '''
    sigma lists
    '''
    sigma_offs_est = []
    sigma_offs_actual = []

    sigma_amp_est = []
    sigma_amp_actual = []

    sigma_peaks_est = []
    sigma_peaks_actual = []

    '''
    freq lists
    '''
    freq_est = []
    freq_actual = []

    '''
    other
    '''
    spikecounts = []
    CV = []
    mu_peaks_est_bad = []

    '''
    Plot the reconstructin of 20 random samples.
    '''
    toplot = arange(numsims)
    shuffle(toplot)
    toplot = toplot[:10]
    for ind, d in enumerate(data.itervalues()):
        i = ind+1
        print "\rProcessing %i of %i" % (i, numsims),
        if not d:
            print "Skipping item %i" % (i)
            continue
        sys.stdout.flush()
        spikes = d.get('spikes')
        mem = d.get('mem')
        mu_offs = d.get('mu_offs')
        mu_amp = d.get('mu_amp')
        sigma_offs = d.get('sigma_offs')
        sigma_amp = d.get('sigma_amp')
        freq = d.get('freq')
        duration = d.get('duration')
        filehead = "%s/m%f-%f_s%f-%f-f%f" % (plot_save_loc,
                mu_offs/(mV/ms),
                mu_amp/(mV/ms),
                sigma_offs/(mV/sqrt(ms)),
                sigma_amp/(mV/sqrt(ms)),
                freq/Hz)

        mu_peak = mu_amp+mu_offs
        sigma_peaks = sigma_amp+sigma_offs
        N = len(spikes)
        for ii in range(N):
            st = spikes.get(ii)
            if len(st) < 10: continue
            '''
            Estimator placeholders
            '''
            mo_est = 0
            ma_est = 0
            mp_est = 0
            so_est = 0
            sa_est = 0
            sp_est = 0
            f_est  = 0
            cv = 0
            V = mem[ii]
            cur_isis = diff(st)

            cv = nt.CV(st)

            '''
            Frequency estimation:
                PSD of spike train
            '''
            f_est = psd_est(st, duration)

            '''
            Estimate whether mu offs is super- or sub-threshold.
            '''
            mr_est = is_super(st, duration, f_est, tau)
            mu_rng_est.append(mr_est)


            '''
            Input mu(t) reconstruction:
                Classical mu estimator on all ISIs (folded on period)
            '''
            nbins = 10
            amte = discrete_mu_recon_avg(st, V, f_est+0.1*f_est, nbins, dt=dt)
            mu_t_est.append(amte)

            '''
            mu_peaks estimation error as a function of freq estimation error
            '''
            mpe_bad = []
            freq_mis_errs = linspace(-0.1, 0.1, 20)
            for fme in freq_mis_errs:
                badf = freq+freq*fme
                badamte = discrete_mu_recon_avg(st, V, badf, nbins, dt=dt)
                mpe_bad.append(max(badamte))


            '''
            Mu peak estimation:
                Average peak mu estimation
            '''
            mp_est = max(amte)
            mo_est = mean(amte)

            ma_est = mp_est-mo_est
            #amtsle = slope_mu_recon(st, V, freq, binwidth, w=2*ms, dt=dt)
            #mu_slope_est.append(amtsle)

            #Vn_mean = mslp**2 * pi/2 # estimation of Vn based on mean
            #sigma_est_i = sqrt(2*Vn_mean / (tau * (1+exp(-2*dt/tau))))

            mu_offs_est.append(mo_est)
            mu_offs_actual.append(mu_offs)

            mu_amp_est.append(ma_est)
            mu_amp_actual.append(mu_amp)

            mu_peaks_est.append(mp_est)
            mu_peaks_actual.append(mu_peak)

            sigma_offs_est.append(so_est)
            sigma_offs_actual.append(sigma_offs)

            sigma_amp_est.append(sa_est)
            sigma_amp_actual.append(sigma_amp)

            sigma_peaks_est.append(sp_est)
            sigma_peaks_actual.append(sigma_peaks)

            freq_est.append(f_est)
            freq_actual.append(freq)

            spikecounts.append(len(st))
            CV.append(cv)
            mu_peaks_est_bad.append(mpe_bad)

            if False and ind in toplot:
                '''
                Plot reconstruction.
                '''
                recon = discrete_mu_recon(st, V, dt=dt)
                period = freq**-1
                figure(num=1, figsize=(8, 6), dpi=100)
                clf()
                T = arange(0*second, 1*second, dt)
                plot(st, recon, linestyle='-', marker='.',
                        color='black', markersize=12, label='$\hat{\mu}$')
                plot(T, sin(T*2*pi*freq)*mu_amp+mu_offs, linestyle='--',
                        color='grey', label='$\mu(t)$')
                axis(xmax=1)
                #legend(fontsize=15)
                legend()
                xlabel('t (sec)', size=15)
                ylabel('mV/ms', size=15)
                xticks(size=15)
                yticks(size=15)
                savefig('%s_mu_t_est_actual.png' % filehead)
                savefig('%s_mu_t_est_actual.pdf' % filehead)

                '''
                Plot `folded` reconstruction.
                '''
                figure(num=1, figsize=(8, 6), dpi=100)
                clf()
                T = arange(0*second, period, dt)
                plot(st % period, recon, linestyle='', marker='.',
                        color='black', markersize=12, label='$\hat{\mu}$')
                plot(T, sin(T*2*pi*freq)*mu_amp+mu_offs, linestyle='--',
                        color='grey', label='$\mu(t)$')
                #legend(fontsize=15)
                legend()
                xlabel('t (sec)', size=15)
                ylabel('mV/ms', size=15)
                xticks(size=15)
                yticks(size=15)
                savefig('%s_mu_t_est_actual_period.png' % filehead)
                savefig('%s_mu_t_est_actual_period.pdf' % filehead)

                '''
                Plot binned reconstruction.
                '''
                figure(num=1, figsize=(8, 6), dpi=100)
                clf()
                T = arange(0*second, period, dt)
                plot(linspace(0*second, period, len(amte)), amte,
                        linestyle='-', marker='.',
                        color='black', markersize=12, label='$\hat{\mu}_b$')
                plot(T, sin(T*2*pi*freq)*mu_amp+mu_offs, linestyle='--',
                        color='grey', label='$\mu(t)$')
                #legend(fontsize=15)
                legend()
                xlabel('t (sec)', size=15)
                ylabel('mV/ms', size=15)
                xticks(size=15)
                yticks(size=15)
                savefig('%s_mu_t_est_actual_binned.png' % filehead)
                savefig('%s_mu_t_est_actual_binned.pdf' % filehead)



    print "\n----\n"
    plot_est_vs_act(mu_peaks_actual, mu_peaks_est,
            '$\mu_p (mV/ms)$', '$\hat{\mu}_p (mV/ms)$', 'mu_p_estimation')
    plot_est_vs_act_err(mu_peaks_actual, mu_peaks_est,
            '$\mu_p (mV/ms)$', '$\hat{\mu}_p (mV/ms)$', 'mu_p_est_err')
    plot_est_vs_act(mu_offs_actual, mu_offs_est,
            '$\mu_0 (mV/ms)$', '$\hat{\mu}_0 (mV/ms)$', 'mu_0_estimation')
    plot_est_vs_act_err(mu_offs_actual, mu_offs_est,
            '$\mu_0 (mV/ms)$', '$\hat{\mu}_0 (mV/ms)$', 'mu_0_est_err')
    plot_est_vs_act(mu_amp_actual, mu_amp_est,
            '$\mu_a (mV/ms)$', '$\hat{\mu}_a (mV/ms)$', 'mu_a_estimation')
    plot_est_vs_act_err(mu_amp_actual, mu_amp_est,
            '$\mu_a (mV/ms)$', '$\hat{\mu}_a (mV/ms)$', 'mu_a_est_err')

    #plot_est_vs_act(sigma_peaks_actual, sigma_peaks_est,
    #        '$\sigma_p$ (mV/sqrt(ms))',
    #        '$\hat{\sigma_p}$ (mV/sqrt(ms))',
    #        'sigma_p_estimation')

    plot_est_vs_act(freq_actual, freq_est, '$f$ (Hz)', '$\hat{f}$ (Hz)',
                                        'frequency_estimation')
    plot_est_vs_act_err(freq_actual, freq_est, '$f$ (Hz)', '$\hat{f}$ (Hz)',
                                        'frequency_est_err')


    fe_err = abs(array(freq_est)-array(freq_actual))/array(freq_actual)
    '''
    CV vs frequency estimation error
    '''
    figure(num=1, figsize=(8, 6), dpi=100)
    clf()
    scatter(CV, fe_err, c='black')
    xlabel('CV', size=15)
    ylabel('$f$ est. error', size=15)
    xticks(size=15)
    yticks(size=15)
    savefig('%s/CV_freq_est_err.png' % plot_save_loc)
    savefig('%s/CV_freq_est_err.pdf' % plot_save_loc)


    '''
    mu offs vs mu actual for non-zero frequency estimation errors
    '''
    figure(num=1, figsize=(8, 6), dip=100)
    clf()
    nzerr = flatnonzero(fe_err)
    scatter(mu_offs_actual, mu_amp_actual, c='black', marker='.', s=10,
            label=r'$ \varepsilon_f  = 0 $')
    scatter(array(mu_offs_actual)[nzerr], array(mu_amp_actual)[nzerr],
            c='black', marker='^', s=100, label=r'$ \varepsilon_f > 0 $')
    xlabel('$\mu_0$ (mV/ms)', size=15)
    ylabel('$\mu_a$ (mV/ms)', size=15)
    xticks(fontsize=15)
    yticks(fontsize=15)
    legend()
    title("Non-zero frequency estimation errors", size=15)
    savefig('%s/mu_for_nz_freq_err.png' % plot_save_loc)
    savefig('%s/mu_for_nz_freq_err.pdf' % plot_save_loc)

    '''
    Up to this point, everything was saved in lists. After the savez command
    they are saved as numpy arrays.
    '''
    archivename = '%s.npz' % data_dir
    np.savez(archivename,
            mu_offs_est=mu_offs_est,
            mu_offs_actual=mu_offs_actual,
            mu_amp_est=mu_amp_est,
            mu_amp_actual=mu_amp_actual,
            mu_peaks_est=mu_peaks_est,
            mu_peaks_actual=mu_peaks_actual,
            mu_t_est=mu_t_est,
            mu_rng_est=mu_rng_est,
            sigma_offs_est=sigma_offs_est,
            sigma_offs_actual=sigma_offs_actual,
            sigma_amp_est=sigma_amp_est,
            sigma_amp_actual=sigma_amp_actual,
            sigma_peaks_est=sigma_peaks_est,
            sigma_peaks_actual=sigma_peaks_actual,
            freq_est=freq_est,
            freq_actual=freq_actual,
            spikecounts=spikecounts,
            CV=CV,
            mu_peaks_est_bad=mu_peaks_est_bad,
            )
    print "Actual and estimated values have been saved in the %s archive" % (
            archivename)


'''
Use the following to load the data in another script or an interactive
session:

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



'''
