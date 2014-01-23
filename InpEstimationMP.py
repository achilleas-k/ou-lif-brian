import matplotlib
matplotlib.use('Agg') # for plotting remotely or in screen (i.e., no DISPLAY)
from matplotlib import rc
from brian import *
from brian.tools.datamanager import *
import sys, os
from time import time, sleep
import neurotools as nt
from multiprocessing import Process, Lock, Queue, Pool
from uuid import uuid4


nproc = 8 # number of concurrent processes
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
    filename_pdf = '%s/%s.pdf' % (plot_save_loc, figname)
    filename_png = '%s/%s.png' % (plot_save_loc, figname)
    savefig(filename_pdf)
    savefig(filename_png)
    print 'Figures saved: %s and %s' % (filename_pdf, filename_png)
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
    rc('text', usetex=True) # use latex
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
    savetypes = ["pdf", "png", "eps", "svg"]
    for ft in savetypes:
        filename = "%s/%s.%s" % (plot_save_loc, figname, ft)
        savefig(filename)
        print "Figures saved: %s" % (filename)
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
    acorr = acorr[len(acorr)/2:] # acorr is symmetric, keep only positive
    psd = fft(acorr)
    #cutfreq = len(psd)/2
    #psd[cutfreq:] = 0 # low pass filter
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

def discrete_sigma_recon(spikes, V, bin=0.1*ms, dt=0.1*ms):
    '''
    Estimates the value of sigma for each interspike interval based on the
    sigma estimator from Bibbona et al., 2008.

    Calls mu estimation again, which is redundant since we calculated it before.
    I should change it so that it's estimated only once.
    '''
    mu_est = discrete_mu_recon(spikes, V, dt=dt)
    recon = []
    st = append(0, spikes)
    for tp, t, mu in zip(st[:-1], st[1:], mu_est):
        t *= second
        tp *= second
        isi = t-tp
        K = round(isi/dt)
        Vi = V[(tp/dt+1):(t/dt+1)]
        se2 = 2/(K-1) * sum(((Vi[1:] - mu*tau + (mu*tau - Vi[:-1]) *\
                exp(-dt/tau))**2) / (tau*(1-exp(-2*dt/tau))))
        se = sqrt(se2)
        recon.append(se)
    return array(recon)

def discrete_sigma_recon_avg(spikes, V, freq, nbins, slopebin=0.1*ms, dt=0.1*ms):
    '''
    Estimates the value of sigma for each interspike interval based on the sigma
    estimator from Bibbona et al., 2008, folds the estimated values
    based on period using the frequency `freq` and averages the estimated
    values in `nbins` bins.

    TODO: Try a moving average instead of binning.
    '''
    slopebin = max(dt, slopebin)
    recon = discrete_sigma_recon(spikes, V, slopebin, dt)
    if is_dimensionless(freq):
        freq *= Hz
    period = 1/freq
    periodspikes = spikes % period
    binedges = linspace(0*second, period, nbins+1)
    avgsigma = []
    for l, r in zip(binedges[:-1], binedges[1:]):
        binspike_inds = bitwise_and(periodspikes > l, periodspikes < r)
        binsigma = recon[binspike_inds]
        avgsigma.append(mean(binsigma) if len(binsigma) else 0)
    return array(avgsigma)

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

def estimateparams(
            d,
            mu_offs_est,
            mu_offs_actual,
            mu_amp_est,
            mu_amp_actual,
            mu_peaks_est,
            mu_peaks_actual,
            mu_t_est,
            mu_rng_est,
            sigma_offs_est,
            sigma_offs_actual,
            sigma_amp_est,
            sigma_amp_actual,
            sigma_peaks_est,
            sigma_peaks_actual,
            freq_actual,
            freq_est,
            spikecounts,
            CV,
            mu_peaks_est_bad,
            i,
            ):

    if not d:
        return

    # used to handle non-trivial unit conversion for sigmas
    sigma_units = mV/sqrt(ms)
    spikes = d.get('spikes')
    mem = d.get('mem')
    mu_offs = d.get('mu_offs')
    mu_amp = d.get('mu_amp')
    sigma_offs = d.get('sigma_offs')/sigma_units
    sigma_amp = d.get('sigma_amp')/sigma_units
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


        '''
        Input mu(t) reconstruction:
            Bibbona mu estimator on all ISIs (folded on period)
        '''
        nbins = 10
        amte = discrete_mu_recon_avg(st, V, f_est, nbins, dt=dt)

        '''
        mu_peaks estimation error as a function of freq estimation error
        '''
        mpe_bad = []
        freq_mis_errs = linspace(-0.1, 0.1, 21)
        for fme in freq_mis_errs:
            badf = freq*(1+fme)
            badamte = discrete_mu_recon_avg(st, V, badf, nbins, dt=dt)
            mpe_bad.append(max(badamte))


        '''
        Mu peak estimation:
            Average peak mu estimation
        Mu offs estimation:
            Mean of amte
        '''
        mp_est = max(amte)
        mo_est = mean(amte)
        ma_est = mp_est-mo_est

        '''
        Input sigma(t) reconstruction:
            Bibbona sigma estimator on all ISIs (folded on period)

        NOTICE: mte was already calculated in discrete_mu_recon_avg but
                in order to avoid changing all the functions for now I recalc
                inside discrete_sigma_recon
        '''
        nbins = 10
        aste = discrete_sigma_recon_avg(st, V, freq, nbins)

        '''
        Sigma peak and offs estimation
        '''
        sp_est = max(aste)
        so_est = max(mean(aste), sp_est/2)
        sa_est = sp_est-so_est

        sp_est /= sigma_units
        so_est /= sigma_units
        sa_est /= sigma_units

        '''Inserting results into queues'''
        mu_offs_est.put((i, mo_est))
        mu_offs_actual.put((i, mu_offs))

        mu_amp_est.put((i, ma_est))
        mu_amp_actual.put((i, mu_amp))

        mu_peaks_est.put((i, mp_est))
        mu_peaks_actual.put((i, mu_peak))

        mu_rng_est.put((i, mr_est))
        mu_t_est.put((i, amte))

        sigma_offs_est.put((i, so_est))
        sigma_offs_actual.put((i, sigma_offs))

        sigma_amp_est.put((i, sa_est))
        sigma_amp_actual.put((i, sigma_amp))

        sigma_peaks_est.put((i, sp_est))
        sigma_peaks_actual.put((i, sigma_peaks))

        freq_est.put((i, f_est))
        freq_actual.put((i, freq))

        spikecounts.put((i, len(st)))
        CV.put((i, cv))
        mu_peaks_est_bad.put((i, mpe_bad))

    return

def dequeue_item(aqueue, anarray):
    if aqueue.empty():
        return
    ind, dat = aqueue.get()
    anarray[ind] = dat
    return

if __name__ == '__main__':
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
    mu_offs_est = [0]*numsims
    mu_offs_actual = [0]*numsims

    mu_amp_est = [0]*numsims
    mu_amp_actual = [0]*numsims

    mu_peaks_est = [0]*numsims
    mu_peaks_actual = [0]*numsims
    mu_t_est = [0]*numsims
    mu_rng_est = [0]*numsims

    '''
    sigma lists
    '''
    sigma_offs_est = [0]*numsims
    sigma_offs_actual = [0]*numsims

    sigma_amp_est = [0]*numsims
    sigma_amp_actual = [0]*numsims

    sigma_peaks_est = [0]*numsims
    sigma_peaks_actual = [0]*numsims

    '''
    freq lists
    '''
    freq_est = [0]*numsims
    freq_actual = [0]*numsims

    '''
    other
    '''
    spikecounts = [0]*numsims
    CV = [0]*numsims
    mu_peaks_est_bad = [0]*numsims

    '''
    mu queues
    '''
    mu_offs_est_q = Queue(numsims)
    mu_offs_actual_q = Queue(numsims)

    mu_amp_est_q = Queue(numsims)
    mu_amp_actual_q = Queue(numsims)

    mu_peaks_est_q = Queue(numsims)
    mu_peaks_actual_q = Queue(numsims)

    mu_t_est_q = Queue(numsims)
    mu_rng_est_q = Queue(numsims)

    '''
    sigma queues
    '''
    sigma_offs_est_q = Queue(numsims)
    sigma_offs_actual_q = Queue(numsims)

    sigma_amp_est_q = Queue(numsims)
    sigma_amp_actual_q = Queue(numsims)

    sigma_peaks_est_q = Queue(numsims)
    sigma_peaks_actual_q = Queue(numsims)

    '''
    freq queues
    '''
    freq_est_q = Queue(numsims)
    freq_actual_q = Queue(numsims)

    '''
    other
    '''
    spikecounts_q = Queue(numsims)
    CV_q = Queue(numsims)
    mu_peaks_est_bad_q = Queue(numsims)

    plist = []
    rlist = []
    nfin = 0
    i = 0
    print "Building process list ..."
    for d in data.itervalues():
        sys.stdout.write('\r%i/%i' % (i+1, numsims))
        sys.stdout.flush()
        plist.append(Process(target=estimateparams,
            kwargs={
                    'd': d,
                    'mu_offs_est': mu_offs_est_q,
                    'mu_offs_actual': mu_offs_actual_q,
                    'mu_amp_est': mu_amp_est_q,
                    'mu_amp_actual': mu_amp_actual_q,
                    'mu_peaks_est': mu_peaks_est_q,
                    'mu_peaks_actual': mu_peaks_actual_q,
                    'mu_t_est': mu_t_est_q,
                    'mu_rng_est': mu_rng_est_q,
                    'sigma_offs_est': sigma_offs_est_q,
                    'sigma_offs_actual': sigma_offs_actual_q,
                    'sigma_amp_est': sigma_amp_est_q,
                    'sigma_amp_actual': sigma_amp_actual_q,
                    'sigma_peaks_est': sigma_peaks_est_q,
                    'sigma_peaks_actual': sigma_peaks_actual_q,
                    'freq_actual': freq_actual_q,
                    'freq_est': freq_est_q,
                    'spikecounts': spikecounts_q,
                    'CV': CV_q,
                    'mu_peaks_est_bad': mu_peaks_est_bad_q,
                    'i': i,
                    }
            ))
        i += 1

    print "\nProcess list created. Starting processes ... "
    while plist or rlist:
        running = 0
        sys.stdout.write('\r%i of %i complete' % (nfin, numsims))
        sys.stdout.flush()
        while len(rlist) < nproc and plist:
            start_proc = plist.pop()
            start_proc.start()
            running += 1
            rlist.append((start_proc, time()))
            #print "Process %i started" % (start_proc.pid)

        for rproc, stime in rlist[:]:
            if not rproc.is_alive() or rproc.exitcode is not None:
                rlist.remove((rproc, stime))
                #print "Process %i finished after %f seconds from start" % (
                #        rproc.pid, time()-stime)
                nfin += 1
                running -= 1
            else:
                running += 1
                #print "Process %i has been running for %f seconds" % (
                #        rproc.pid, time()-stime)

        #print '-------------------'
        'Dequeue here'
        dequeue_item(mu_peaks_est_q, mu_peaks_est)
        dequeue_item(mu_peaks_actual_q, mu_peaks_actual)
        dequeue_item(mu_offs_est_q, mu_offs_est)
        dequeue_item(mu_offs_actual_q, mu_offs_actual)
        dequeue_item(mu_amp_est_q, mu_amp_est)
        dequeue_item(mu_amp_actual_q, mu_amp_actual)
        dequeue_item(mu_t_est_q, mu_t_est)
        dequeue_item(mu_rng_est_q, mu_rng_est)
        dequeue_item(sigma_peaks_est_q, sigma_peaks_est)
        dequeue_item(sigma_peaks_actual_q, sigma_peaks_actual)
        dequeue_item(sigma_offs_est_q, sigma_offs_est)
        dequeue_item(sigma_offs_actual_q, sigma_offs_actual)
        dequeue_item(sigma_amp_est_q, sigma_amp_est)
        dequeue_item(sigma_amp_actual_q, sigma_amp_actual)
        dequeue_item(freq_est_q, freq_est)
        dequeue_item(freq_actual_q, freq_actual)
        dequeue_item(spikecounts_q, spikecounts)
        dequeue_item(CV_q, CV)
        dequeue_item(mu_peaks_est_bad_q, mu_peaks_est_bad)

    sys.stdout.write('\n------------------\n')

    print "Cleaning arrays..."
    nzind = flatnonzero(array(mu_peaks_actual))

    mu_peaks_est = [mu_peaks_est[i] for i in nzind]
    mu_peaks_actual = [mu_peaks_actual[i] for i in nzind]
    mu_offs_est = [mu_offs_est[i] for i in nzind]
    mu_offs_actual = [mu_offs_actual[i] for i in nzind]
    mu_amp_est = [mu_amp_est[i] for i in nzind]
    mu_amp_actual = [mu_amp_actual[i] for i in nzind]
    mu_t_est = [mu_t_est[i] for i in nzind]
    mu_rng_est = [mu_rng_est[i] for i in nzind]

    sigma_peaks_est = [sigma_peaks_est[i] for i in nzind]
    sigma_peaks_actual = [sigma_peaks_actual[i] for i in nzind]
    sigma_offs_est = [sigma_offs_est[i] for i in nzind]
    sigma_offs_actual = [sigma_offs_actual[i] for i in nzind]
    sigma_amp_est = [sigma_amp_est[i] for i in nzind]
    sigma_amp_actual = [sigma_amp_actual[i] for i in nzind]

    freq_est = [freq_est[i] for i in nzind]
    freq_actual = [freq_actual[i] for i in nzind]
    CV = [CV[i] for i in nzind]
    mu_peaks_est_bad = [mu_peaks_est_bad[i] for i in nzind]

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


    print "\n----\n"
    print "Generating figures ..."
    plot_est_vs_act(freq_actual, freq_est, '$f$ (Hz)', '$\hat{f}$ (Hz)',
                                        'frequency_estimation', plot_zeros=True)
    plot_est_vs_act_err(freq_actual, freq_est, '$f$ (Hz)', '$\hat{f}$ (Hz)',
                                        'frequency_est_err', plot_zeros=True)

    plot_est_vs_act(mu_peaks_actual, mu_peaks_est,
            '$\mu_p (mV/ms)$', '$\hat{\mu}_p (mV/ms)$', 'mu_p_estimation',
            plot_zeros=True)
    plot_est_vs_act_err(mu_peaks_actual, mu_peaks_est,
            '$\mu_p (mV/ms)$', '$\hat{\mu}_p (mV/ms)$', 'mu_p_est_err',
            plot_zeros=True)
    plot_est_vs_act(mu_offs_actual, mu_offs_est,
            '$\mu_0 (mV/ms)$', '$\hat{\mu}_0 (mV/ms)$', 'mu_0_estimation',
            plot_zeros=True)
    plot_est_vs_act_err(mu_offs_actual, mu_offs_est,
            '$\mu_0 (mV/ms)$', '$\hat{\mu}_0 (mV/ms)$', 'mu_0_est_err',
            plot_zeros=True)
    plot_est_vs_act(mu_amp_actual, mu_amp_est,
            '$\mu_a (mV/ms)$', '$\hat{\mu}_a (mV/ms)$', 'mu_a_estimation',
            plot_zeros=True)
    plot_est_vs_act_err(mu_amp_actual, mu_amp_est,
            '$\mu_a (mV/ms)$', '$\hat{\mu}_a (mV/ms)$', 'mu_a_est_err',
            plot_zeros=True)

    plot_est_vs_act(sigma_peaks_actual, sigma_peaks_est,
            '$\sigma_p (mV/\sqrt{ms})$',
            '$\hat{\sigma}_p (mV/\sqrt{ms})$',
            'sigma_p_estimation', plot_zeros=True)
    plot_est_vs_act_err(sigma_peaks_actual, sigma_peaks_est,
            '$\sigma_p (mV/\sqrt{ms})$',
            '$\hat{\sigma}_p (mV/\sqrt{ms})$',
            'sigma_p_est_err', plot_zeros=True)
    plot_est_vs_act(sigma_offs_actual, sigma_offs_est,
            '$\sigma_0 (mV/\sqrt{ms})$',
            '$\hat{\sigma}_0 (mV/\sqrt{ms})$',
            'sigma_0_estimation', plot_zeros=True)
    plot_est_vs_act_err(sigma_offs_actual, sigma_offs_est,
            '$\sigma_0 (mV/\sqrt{ms})$',
            '$\hat{\sigma}_0 (mV/\sqrt{ms})$',
            'sigma_0_est_err', plot_zeros=True)
    plot_est_vs_act(sigma_amp_actual, sigma_amp_est,
            '$\sigma_a (mV/\sqrt{ms})$',
            '$\hat{\sigma}_a (mV/\sqrt{ms})$',
            'sigma_a_estimation', plot_zeros=True)
    plot_est_vs_act_err(sigma_amp_actual, sigma_amp_est,
            '$\sigma_a (mV/\sqrt{ms})$',
            '$\hat{\sigma}_a (mV/\sqrt{ms})$',
            'sigma_a_est_err', plot_zeros=True)


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
    mu offs vs mu amp for non-zero frequency estimation errors
    '''
    figure(num=1, figsize=(8, 6), dpi=100)
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
    mu peaks error vs freq misestimation
    '''
    figure(num=1, figsize=(8, 6), dpi=100)
    clf()
    bad_errs = array([abs(mpeb-mu_peaks_actual)/mu_peaks_actual
                                for mpeb in transpose(mu_peaks_est_bad)])
    avg_bad_errs = mean(bad_errs, axis=1)
    fme = linspace(-0.1, 0.1, len(avg_bad_errs))
    plot(fme, avg_bad_errs, color='black')
    axis(xmin=min(fme), xmax=max(fme))
    xlabel(r'$\alpha_f$', size=15)
    ylabel(r'$\varepsilon_{\mu_p}$', size=15)
    yticks(fontsize=15)
    xticks(fontsize=15)
    savefig('%s/mu_err_vs_freq_err.png' % plot_save_loc)
    savefig('%s/mu_err_vs_freq_err.pdf' % plot_save_loc)


    print "ALL DONE!"

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
