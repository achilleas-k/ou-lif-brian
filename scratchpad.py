import sys
print ("This is a scratchpad. It contains code segments that are not meant to "
"be run together (or at all). It is meant for copying and pasting code into "
"other scripts or more commonly the python shell. "
"\nExiting!")
sys.exit(0)


from scipy.stats import norm, kstest
v_th = 10*mV
tau = 10*ms

def sinewave(offs, amp, freq, t):
    return offs+amp*sin(2*pi*freq*t)

'''
Ditlevsen estimators
'''
def ditlevsen_est(data):
    for d in data.itervalues():
        spikes = d.get('spikes')
        if len(spikes) > 10:
            mu = d.get('mu_offs')
            sigma = d.get('sigma_offs')
            if sigma > 4*mV/sqrt(ms) and mu*tau < v_th:
                theta = (v_th-mu*tau) / (sigma*sqrt(tau))
                th_est = est_theta(spikes)[0]
                thetas = append(thetas, [[theta, th_est]], axis=0)
                thvar = th_est**2 / (len(spikes)*(1-2*th_est**2)**2)
                thvars = append(thvars, thvar)
            if mu*tau > v_th-v_th_eps and mu*tau < v_th+v_th_eps:
                sigma_est = sqrt(mean(2*v_th**2 / (
                    tau*(exp(2*diff(spikes)/tau)-1))))
                sigma_ests = append(sigma_ests, [[sigma, sigma_est]], axis=0)
                svar = sigma_est**2 / (2*len(spikes))
                svars = append(svars, svar)
        else:
            print "Too few spikes %i" % len(spikes)


'''
Firing mu and sigma calculations and plotting (no estimation)
'''
def plot_firing_ms(data):
    for d in data.values():
        mem = d.get('mem')
        spikes = d.get('spikes')
        slopes = d.get('slopes')
        mu_offs = d.get('mu_offs')*ms/mV
        mu_amp = d.get('mu_amp')*ms/mV
        sigma_offs = d.get('sigma_offs')*sqrt(ms)/mV
        sigma_amp = d.get('sigma_amp')*sqrt(ms)/mV
        freq = d.get('freq')/Hz
        firing_mu = []
        firing_sigma = []
        for st in spikes.values():
            firing_mu.extend(sinewave(mu_offs, mu_amp, freq, st))
            firing_sigma.extend(sinewave(sigma_offs, sigma_amp, freq, st))

        filehead = 'tmp/%f_%f_%f_%f_%f' % (
            mu_offs,
            mu_amp,
            sigma_offs,
            sigma_amp,
            freq)

        figure(1)
        h,b,p = hist(firing_mu, bins=50, normed=True)
        mean_mu = mean(firing_mu)
        plot([mean_mu, mean_mu], [0, max(h)], color='red', linewidth=2)
        savefig('%s_firing_mu.png' % (filehead))
        clf()
        figure(2)
        h,b,p = hist(firing_sigma, bins=50, normed=True)
        mean_sigma = mean(firing_sigma)
        plot([mean_sigma, mean_sigma], [0, max(h)], color='red', linewidth=2)
        savefig('%s_firing_sigma.png' % (filehead))
        clf()

'''
Spinner for loading or processing.
'''
def spinner():
    i=0
    chars='/-\\|'
    lc = len(chars)
    while 1:
        i+=1
        yield chars[i % lc]

'''
Calculate maximum firing rate based on moving average of binned spike train
'''
import neurotools as nt
def max_rate_movavg(data):
    maxfrate = array([])
    meanfrate = array([])
    mu_as = array([])
    mu_os = array([])
    for d in data.itervalues():
        freq = d.get('freq')
        mu_offs = d.get('mu_offs')
        mu_amp = d.get('mu_amp')
        for st in spikes.itervalues():
            bst = nt.times_to_bin(st)
            frate = movavg(bst, 1/(freq*dt))
            maxfrate = append(maxfrate, max(frate))
            mu_as = append(mu_as, mu_amp)
            mu_os = append(mu_os, mu_offs)
            print ".",

'''
Ditlevsen estimator sub threshold (Z parameter)
'''
def est_mu_Z(data):
    mu_offs, mu_amp, frate_peak, frate_fmean, frate_mean = [], [], [], [], []
    Z = []
    tau = 10*ms
    v_th = 10*mV
    for d in data.itervalues():
        isis = []
        spikes = d.get('spikes')
        freq = d.get('freq')
        for st in spikes.itervalues():
            if len(st) < 2:
                continue
            isis.extend(diff(st))
        if not len(isis):
            continue
        h,b = histogram(isis, bins=100, density=True)
        halfperiod = 0.5/freq
        filt_isis = [fisi for fisi in isis if fisi < float(halfperiod)]
        frate_peak.append(b[h == max(h)][0])
        frate_fmean.append(mean(filt_isis))
        frate_mean.append(mean(isis))
        mu_amp.append(d.get('mu_amp'))
        mu_offs.append(d.get('mu_offs'))
        Z.append(mean(exp(array(isis)/tau)))
        mu_est = v_th*array(Z)/(tau*(array(Z)-1))


'''
I think this is the same as above (???)
'''
def est_mu_Z2(data):
    mu_offs, mu_amp, mean_peak_rate, mu_est = [],[],[],[]
    for d in data.itervalues():
        spikes = d.get('spikes')
        freq = d.get('freq')
        if len(spikes) == 0:
            continue
        mrap, Z = mean_rate_at_peak(spikes, freq)
        if mrap:
            mean_peak_rate.append(mrap)
        else:
            continue
        mu_offs.append(d.get('mu_offs'))
        mu_amp.append(d.get('mu_amp'))
        Zarr = array(Z)
        tau = 10*ms
        mu_est.append((v_th*Zarr)/(tau*(Zarr-1)))


'''
Estimates frequency of input sine wave, plots the power spectrum marked with
estimated and actual frequencies and the spike train (as a stem plot) with
superimposed estimated and actual input sine waves.
'''
from brian import *; from brian.tools.datamanager import *
def est_plot_freq(data):
    for d in data.itervalues():
        spikes = d.get('spikes')
        mu_offs = d.get('mu_offs')
        mu_amp = d.get('mu_amp')
        sigma_offs = d.get('sigma_offs')
        sigma_amp = d.get('sigma_amp')
        freq_actual = d.get('freq')
        print '==================='
        for i, st in enumerate(spikes.values()):
            if len(st) < 2:
                continue
            bin = 10*ms
            dt = 0.1*ms
            width = st[-1]*0.75
            acorr = autocorrelogram(st, width=width, bin=bin)
            acorr = acorr[len(acorr)/2:]
            powersp = fft(acorr)
            cutfreq = len(powersp)/2
            powersp[cutfreq:] = zeros(len(powersp[cutfreq:]))
            powersp0 = append(0, powersp[1:])
            try:
                max_amp_i = int(min(argmax(powersp0)))
            except TypeError:
                max_indices = argmax(powersp0)
                if (len(shape(max_indices))) > 1:
                        max_indices = reshape(max_indices,
                                        (product(shape(max_indices))))
                max_amp_i = min(max_indices)
            max_freq = float(max_amp_i)/(len(powersp0)*bin)
            freq_est = max_freq
            relerr = abs(freq_actual-freq_est)/freq_actual

            t = frange(0, 3, 0.001)
            actual_wave = 1+0.5*sin(t*2*pi*freq_actual)
            est_wave = 1+0.5*sin(t*2*pi*freq_est)
            filehead = "tmp/act%f-est%f-err%f-%i" % (
                        freq_actual,
                        freq_est,
                        relerr,
                        i)
            plot(powersp, color='blue', label='power')
            plot([max_amp_i, max_amp_i], [0, powersp[max_amp_i]],
                    color='red',
                    linestyle='--',
                    marker='o',
                    linewidth=2,
                    label='estimate')
            actual_i = int(freq_actual*len(powersp)*bin)
            plot([actual_i, actual_i], [0, powersp[actual_i]],
                    color='green',
                    linestyle='--',
                    marker='o',
                    linewidth=2,
                    label='actual')
            axis(ymax=max(powersp0)*1.1, ymin=min(powersp0)*1.1)
            legend(loc='best')
            title('Actual freq: %s, Estimate: %s\n'
                    'Mo: %s, Ma: %s\n'
                    'So: %s, Sa: %s\n' % (
                        freq_actual, freq_est,
                        display_in_unit(mu_offs, mV/ms),
                        display_in_unit(mu_amp, mV/ms),
                        display_in_unit(sigma_offs, mV/sqrt(ms)),
                        display_in_unit(sigma_amp, mV/sqrt(ms))))
            subplots_adjust(top=0.85)
            filename = '%s-power.png' % (filehead)
            savefig(filename)
            clf()
            stem(st, ones(len(st)), markerfmt='.', linefmt='b-')
            plot(t, actual_wave, color='green', label='Actual')
            plot(t, est_wave, color='red', label='Estimated')
            legend(loc='best')
            axis(ymax=2)
            title('Actual freq: %s, Estimate: %s\n'
                    'Mo: %s, Ma: %s\n'
                    'So: %s, Sa: %s\n' % (
                        freq_actual, freq_est,
                        display_in_unit(mu_offs, mV/ms),
                        display_in_unit(mu_amp, mV/ms),
                        display_in_unit(sigma_offs, mV/sqrt(ms)),
                        display_in_unit(sigma_amp, mV/sqrt(ms))))
            subplots_adjust(top=0.85)
            filename = '%s-spikes.png' % (filehead)
            savefig(filename)
            clf()
            print "A:%s, E:%s" % (freq_actual, freq_est)

'''
Plot spike train as stem plot and superimpose moving average.
See which is the best averaging window.
'''
from neurotools import times_to_bin
def plot_spikes_movavg(data):
    dt = float(0.1*ms)
    t = frange(0, 3, dt)
    for d in data.itervalues():
        spikes = d.get('spikes')
        mu_offs = d.get('mu_offs')
        mu_amp = d.get('mu_amp')
        sigma_offs = d.get('sigma_offs')
        sigma_amp = d.get('sigma_amp')
        freq = d.get('freq')
        filehead = "tmp/m%f-%f_s%f-%f-f%f" % (mu_offs/(mV/ms),
                        mu_amp/(mV/ms),
                        sigma_offs/(mV/sqrt(ms)),
                        sigma_amp/(mV/sqrt(ms)),
                        freq/Hz)

        for i, st in enumerate(spikes.values()):
            if len(st) < 5:
                continue
            binst = times_to_bin(st)
            avgwin = 1.5*mean(diff(st))/dt
            movrate = movavg(binst, avgwin)
            stem(st, ones(len(st)), markerfmt='')
            plot(t[:len(movrate)]+0.5*avgwin*dt, movrate/mean(movrate)*1.1,
                    color='red', label='Rate')
            plot(t, 1+0.5*sin(t*2*pi*freq), color='green', label='Inp. wave')
            filename = "%s-%i-avgrate.png" % (filehead, i)
            axis(xmax=1)
            savefig(filename)
            clf()

'''
Find best averaging window by taking moving averages of the spike train using
a range of values and using them to estimate the underlying frequency.
Minimising the estimation error and observing the best window values should
provide some insight.
'''
from neurotools import times_to_bin
import sys
import os
def movavg_rate_error(data):
    dt = float(0.1*ms)
    t = frange(0, 3, dt)
    freq_act = []
    freq_est = []
    best_wins = []
    for d in data.itervalues():
        spikes = d.get('spikes')
        mu_offs = d.get('mu_offs')
        mu_amp = d.get('mu_amp')
        sigma_offs = d.get('sigma_offs')
        sigma_amp = d.get('sigma_amp')
        freq = d.get('freq')
        filehead = "tmp/m%f-%f_s%f-%f-f%f" % (mu_offs/(mV/ms),
                        mu_amp/(mV/ms),
                        sigma_offs/(mV/sqrt(ms)),
                        sigma_amp/(mV/sqrt(ms)),
                        freq/Hz)
        print("="*100)

        for i, st in enumerate(spikes.values()):
            if len(st) < 5:
                continue
            binst = times_to_bin(st)
            isis = diff(st)
            winscales = linspace(1*ms, 100*ms, 100)/dt
            winscales.astype('int')
            all_win_estimates = []
            for ws in winscales:
                sys.stdout.write('.')
                sys.stdout.flush()
                x = frange(-ws*4, ws*4, 1) # 4 stds per side
                win = normpdf(x, 0, ws)
                frate = convolve(binst, win, mode='valid')
                fft_rate = fft(frate)
                fft0 = append(0, fft_rate[1:])
                max_amp_i = argmax(fft0)
                freq_est_cur = max_amp_i/(len(fft0)*dt)
                all_win_estimates.append(freq_est_cur)

            all_win_errors = abs(subtract(all_win_estimates,freq))
            best_est_i = argmin(all_win_errors)
            best_est = all_win_estimates[best_est_i]
            best_win = winscales[best_est_i]
            freq_est.append(best_est)
            freq_act.append(freq)
            best_wins.append(best_win)
            sys.stdout.write('\rAct: %f, Est: %f [%f, %f], Win: %f [%f, %f]\n' % (
                freq, best_est, min(all_win_estimates), max(all_win_estimates),
                best_win, min(winscales), max(winscales)))

'''
Plot PSDs of all the spike trains which resulted in an erroneous estimation.
PSDs which include NaN or inf are marked.
'''
def plot_bad_PSD(err_ind, spiketrains):
    ioff()
    for ii in err_ind:
        st = spiketrains[ii]
        width = 1
        bin = 1*ms
        acorr = autocorrelogram(st, width=width, bin=bin)
        acorr = acorr[len(acorr)/2:]
        fft_acorr = fft(acorr)
        cutfreq = len(fft_acorr)/2
        fft_acorr[cutfreq:] = 0
        fft_acorr[0] = 0
        est_freq_i = argmax(fft_acorr)
        est_freq = (1.0*est_freq_i/len(fft_acorr))/bin
        act_freq = freq[ii]*Hz
        act_freq_i = int(act_freq*len(fft_acorr)*bin)
        title_str = ('est freq: %s, act freq: %s' % (est_freq, act_freq))
        if ii in nan_ind:
            title_str += ' inf or nan in PSD'
            line_height = 1
        else:
            line_height = max(fft_acorr)+1
        plot(fft_acorr, label='PSD')
        plot([est_freq_i]*2, [-1, line_height], color='red', label='est')
        plot([act_freq_i]*2, [-1, line_height], color='green', label='act')
        legend(loc='best')
        title(title_str)
        axis(xmin=-100)
        savefig('%i.png' % (ii))
        clf()
        print ii
    ion()

'''
Plot spiketrains with the same underlying frequency in a raster, with icreasing
noise level
'''
def plot_strains_by_freq(spiketrains, targetfr, allfr):
    ioff()
    fr = targetfr
    freq = allfr
    common_fr_ind = flatnonzero(freq == fr)
    figure()
    plotted_so = []
    for ii in common_fr_ind:
        so = sigma_offs[ii]*sqrt(ms)/mV
        if so not in plotted_so:
            plotted_so.append(so)
        else:
            continue
        st = spiketrains[ii]
        plot(st, ones(len(st))*so, '.', color='blue')
    ion()



def sinsegments(n):
    clf()
    N = linspace(0, 1, n)
    peakpos = N[n/4]
    t = linspace(0, 1, 1000)
    plot(t, zeros(len(t)), color='black')
    plot(t, sin(t*2*pi), color='blue')
    plot(N, sin(N*2*pi), 'r*')
    plot(peakpos, sin(peakpos*2*pi), 'r*', markersize=25)
    title("N = %i" % n)

    dN = mean(diff(N))
    bin_edges = append(0, N[:-1]+dN/2)
    for be in bin_edges:
        plot([be, be], [-1, 1], 'r--')


def doscatter(filename, ax):
    archive = np.load(filename)
    mu_amp = archive['mu_amp_actual']
    phase_spikes = archive['phase_spikes']
    mu_offs = archive['mu_offs_actual'][0]
    xx, yy = nonzero(phase_spikes)
    xx = mu_amp[xx]
    scat = ax.scatter3D(xx, yy, zs=ones(shape(xx))*mu_offs)
    ax.axis([0,2,0,100])
    return scat


def dowhatev(filename):
    archive = np.load(filename)
    mu_amp = archive['mu_amp_actual']
    mu_offs = archive['mu_offs_actual']
    spikecounts = archive['spikecounts']
    corr = corrcoef(mu_offs, spikecounts)
    scatter(spikecounts, mu_offs)
    title('mu offs: %f, corr: %f' % (mu_offs[0], corr[0,1]))


'''
Category identification
'''
def identify_mo_cat(phase_spikes, mu_amp_actual, mu_offs_actual):
    cat_actual = []
    cat_est = []
    for ps, ma, mo in zip(phase_spikes, mu_amp_actual, mu_offs_actual):
        ps_len = len(ps)
        peak_half = ps[0:ps_len*0.5]
        low_half = ps[ps_len*0.5:]
        min_point = ps[ps_len*0.80:ps_len*0.70]
        first_spike_t = flatnonzero(ps)[0]
        cur_est = 0
        if sum(low_half) == 0:
            cur_est = 1
        elif sum(low_half) > 0:
            cur_est = 3
        cat_est.append(cur_est)


        mp = ma+mo
        mm = mo-ma
        if mo <= 1:
            cat_actual.append(1)
        else:
            cat_actual.append(3)

    cat_est = array(cat_est)
    cat_actual = array(cat_actual)
    errors = abs(cat_est-cat_actual)/cat_actual


'''
Plot categories
'''
def plot_mo_cat(phase_spikes, mu_amp_actual, mu_offs_actual):
    subth = []
    superth = []
    for ps, ma, mo in zip(phase_spikes, mu_amp_actual, mu_offs_actual):
        mp = ma+mo
        mm = mo-ma
        if mo <= 1:
            subth.append(ps)
        else:
            superth.append(ps)

    figure(1)
    title('Sub threshold')
    imshow(subth, aspect='auto', interpolation='none', origin='bottom')
    colorbar()
    figure(2)
    title('Super threshold')
    imshow(superth, aspect='auto', interpolation='none', origin='bottom')
    colorbar()

'''
Calculate binned spike counts after folding
'''
def binned_spike_counts(d):
    freq = d.get('freq')
    period = 1/freq
    duration = d.get('duration')
    nperiods = duration/period
    spikes = d.get('spikes')[0]
    nsegments = 10
    binedges = frange(0*second, duration, period/nsegments)
    binspikes = []
    for l,r in zip(binedges[:-1], binedges[1:]):
        binspikes.append(count_nonzero(spikes[bitwise_and(spikes > l,
            spikes < r)]))
    return array(binspikes)

'''
Estimates the length of each spike burst based on the mean ISI.
A burst is estimated as any series of spikes with intervals < mean ISI.
'''
def find_burst_length(spiketrain):
    '''
    Find the burst length, i.e., the mean consecutive number of spikes with
    an interval < mean ISI.
    '''
    if len(spiketrain) < 2:
        return 0
    ISIs = diff(spiketrain)
    mISI = mean(ISIs)
    sub_mISI = flatnonzero(ISIs < mISI)
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

'''
Reconstruct the input based on the spikecounts in each part of the period.
'''
def reconstruct_inp(d):
    freq = d.get('freq')
    period = 1/freq
    duration = d.get('duration')
    nperiods = duration/period
    spikes = d.get('spikes')[0]
    nsegments = 4
    binwidth = period/nsegments
    binedges = frange(0*second, duration, binwidth)
    recon = []
    for l,r in zip(binedges[:-1], binedges[1:]):
        binspikes = spikes[bitwise_and(spikes > l, spikes < r)]
        spikecount = len(binspikes)
        #if spikecount > 1:
        #    spikewidth = binspikes[-1]-binspikes[0]
        #    recon.append(spikecount/spikewidth)
        #else:
        #    recon.append(0)
        recon.append(spikecount/binwidth)
    return array(recon)/100

'''
Estimate mu based on the spike count reconstruction.
'''
def est_mu(d):
    freq = d.get('freq')
    period = 1/freq
    duration = d.get('duration')
    nperiods = duration/period
    spikes = d.get('spikes')[0]
    nsegments = 2
    binwidth = period/nsegments
    binedges = frange(0*second, duration, binwidth)
    recon = []
    for l,r in zip(binedges[:-1], binedges[1:]):
        binspikes = spikes[bitwise_and(spikes > l, spikes < r)]
        spikecount = len(binspikes)
        if spikecount > 1:
            spikewidth = binspikes[-1]-binspikes[0]
            recon.append(spikecount/spikewidth)
        else:
            recon.append(0)
    recon = array(recon)/100
    mu_peak = max(recon)
    mu_offs = mean(recon)
    mu_amp = mu_peak - mu_offs
    return mu_offs, mu_amp

'''
Calculate the mean interspike interval of all the spikes in a particular
segment of the period, defined as n/N.
'''
def get_misi_segments(d, n, N):
    if n < 1 or type(n) != int:
        raise TypeError('Segment number must be positive integer.')
    if n > N:
        raise Exception('Segment number greater than total segments.')
    freq = d.get('freq')
    preriod = 1/freq
    duration = d.get('duration')
    spikes = d.get('spikes')[0]
    period_edges = frange(0*second, duration, period)
    binned_isis = []
    for pstart in period_edges:
        l = pstart*second+(n-1)/N*period
        r = pstart*second+n/N*period
        binspikes = spikes[bitwise_and(spikes > l, spikes < r)]
        if len(binspikes):
            binned_isis.extend(diff(binspikes))
    return mean(binned_isis)


'''
FULL PROCEDURE SHOULD BE FINALISED HERE
'''
def full_procedure(data):
    vth = 10*mV
    tau = 10*ms
    dt = 0.1*ms
    mp_act = []
    mp_est = []
    mo_est_rng = []
    mo_est = []
    mo_act = []
    mw_ratios = []
    ma_act = []
    for d in data.itervalues():
        mo = d.get('mu_offs')
        ma = d.get('mu_amp')
        mp = mo+ma
        freq = d.get('freq')
        period = 1/freq
        duration = d.get('duration')
        spikes = d.get('spikes')[0]
        V = d.get('mem')[0]

        estimator = 2  # range estimator chooser
        if estimator == 0:
            # mu offs > or < vth
            periodspikes = spikes % period
            highend = 0.35*period
            lowstart = 0.65*period
            nhighspikes = count_nonzero(periodspikes < highend)
            midspikes = count_nonzero(bitwise_and(highend < periodspikes,
                        periodspikes < lowstart))
            lowspikes = count_nonzero(periodspikes > lowstart)
            nperiods = duration/period
            if lowspikes:
                mo_est_rng.append(1)
            elif midspikes:
                firstspike = sorted(periodspikes)[0]*second
                if firstspike > 10*ms: # threshold should be based on vth/tau
                    mo_est_rng.append(0)
                else:
                    mo_est_rng.append(1)
            else:
                mo_est_rng.append(0)
            mo_act.append(mo)
        elif estimator == 1:
            # alt estimator
            if len(spikes):
                periodspikes = spikes % period
                srt_ps = sorted(periodspikes)
                active_phase_width = (srt_ps[-1]-srt_ps[0])*second
                active_phase_ratio = active_phase_width/period
                if active_phase_ratio > 0.35:
                    mo_est_rng.append(1)
                else:
                    mo_est_rng.append(0)
            else:
                mo_est_rng.append(0)
            mo_act.append(mo)
        elif estimator == 2:
            # yet another estimator
            if len(spikes):
                nperiods = int(duration/period)
                widths = []
                for p in range(nperiods):
                    l = p*period
                    r = (p+1)*period
                    hr = r/2 #half period
                    this_period_spikes = spikes[bitwise_and(spikes > l,
                                                                    spikes < hr)]
                    if len(this_period_spikes):
                        this_width = this_period_spikes[-1]-this_period_spikes[0]
                    else:
                        this_width = 0
                    widths.append(this_width)

                mwidth = mean(widths)*second
                mwidth_ratio = mwidth/period
                mw_ratios.append(mwidth_ratio)
                if mwidth_ratio > 0.35:
                    mo_est_rng.append(1)
                else:
                    mo_est_rng.append(0)
            else:
                mw_ratios.append(0)
                mo_est_rng.append(0)
            mo_act.append(mo)

        # estimating mu peak
        if len(spikes) > 2:
            isi = diff(spikes)
            minisi = min(isi)
            minisi_index = argmin(isi)
            min_iii = range((spikes[minisi_index]/dt).astype(int),
                    (spikes[minisi_index+1]/dt).astype(int))
            min_iii = array(min_iii)+1
            K_min = minisi/dt
            mpe = vth/(tau*K_min*(1-exp(-dt/tau))) +\
                                            1/(tau*K_min)*sum(V[min_iii])*volt
            mp_est.append(mpe)

            # mu offs estimation as part of mu peak est
            if mo_est_rng[-1]:
                base_isi = get_misi_segments(d, 1, 4)
          # base_isi_index = argmin(abs(base_isi-isi))
          # base_iii = range((spikes[base_isi_index]/dt).astype(int),
          #         (spikes[base_isi_index+1]/dt).astype(int))
          # base_iii = array(base_iii)+1
          # K_base = base_isi/dt
          # moe = vth/(tau*K_base*(1-exp(-dt/tau))) +\
          #                                 1/(tau*K_base)*sum(V[base_iii])*volt
                mo_est.append(vth/base_isi)
            else:
                mo_est.append(0)
        else:
            mp_est.append(0)
            mo_est.append(0)
        mp_act.append(mp)
        ma_act.append(ma)

    mp_act = array(mp_act)
    mp_est = array(mp_est)
    mo_est_rng = array(mo_est_rng)
    mo_est = array(mo_est)
    mo_act = array(mo_act)
    mw_ratios = array(mw_ratios)

def get_misi_segments(d, n, N):
    if n < 1 or type(n) != int:
        raise TypeError('Segment number must be positive integer.')
    if n > N:
        raise Exception('Segment number greater than total segments.')
    freq = d.get('freq')
    preriod = 1/freq
    duration = d.get('duration')
    spikes = d.get('spikes')[0]
    period_edges = frange(0*second, duration, period)
    binned_isis = []
    for pstart in period_edges:
        l = pstart*second+(n-1)/N*period
        r = pstart*second+n/N*period
        binspikes = spikes[bitwise_and(spikes > l, spikes < r)]
        if len(binspikes):
            binned_isis.extend(diff(binspikes))
    return mean(binned_isis)


'''
FULL PROCEDURE SHOULD BE FINALISED HERE
I have no idea how this differs from above (run a diff)
'''
def full_procedure_alt(data):
    vth = 10*mV
    tau = 10*ms
    dt = 0.1*ms
    mp_act = []
    mp_est = []
    mo_est_rng = []
    mo_est = []
    mo_act = []
    mw_ratios = []
    for d in data.itervalues():
        mo = d.get('mu_offs')
        ma = d.get('mu_amp')
        mp = mo+ma
        freq = d.get('freq')
        period = 1/freq
        duration = d.get('duration')
        spikes = d.get('spikes')[0]
        V = d.get('mem')[0]

        estimator = 2  # range estimator chooser
        if estimator == 0:
            # mu offs > or < vth
            periodspikes = spikes % period
            highend = 0.35*period
            lowstart = 0.65*period
            nhighspikes = count_nonzero(periodspikes < highend)
            midspikes = count_nonzero(bitwise_and(highend < periodspikes,
                        periodspikes < lowstart))
            lowspikes = count_nonzero(periodspikes > lowstart)
            nperiods = duration/period
            if lowspikes:
                mo_est_rng.append(1)
            elif midspikes:
                firstspike = sorted(periodspikes)[0]*second
                if firstspike > 10*ms: # threshold should be based on vth/tau
                    mo_est_rng.append(0)
                else:
                    mo_est_rng.append(1)
            else:
                mo_est_rng.append(0)
            mo_act.append(mo)
        elif estimator == 1:
            # alt estimator
            if len(spikes):
                periodspikes = spikes % period
                srt_ps = sorted(periodspikes)
                active_phase_width = (srt_ps[-1]-srt_ps[0])*second
                active_phase_ratio = active_phase_width/period
                if active_phase_ratio > 0.35:
                    mo_est_rng.append(1)
                else:
                    mo_est_rng.append(0)
            else:
                mo_est_rng.append(0)
            mo_act.append(mo)
        elif estimator == 2:
            # yet another estimator
            if len(spikes):
                nperiods = int(duration/period)
                widths = []
                for p in range(nperiods):
                    l = p*period
                    r = (p+1)*period
                    hr = r/2 #half period
                    this_period_spikes = spikes[bitwise_and(spikes > l,
                                                                    spikes < hr)]
                    if len(this_period_spikes):
                        this_width = this_period_spikes[-1]-this_period_spikes[0]
                    else:
                        this_width = 0
                    widths.append(this_width)

                mwidth = mean(widths)*second
                mwidth_ratio = mwidth/period
                mw_ratios.append(mwidth_ratio)
                if mwidth_ratio > 0.35:
                    mo_est_rng.append(1)
                else:
                    mo_est_rng.append(0)
            else:
                mw_ratios.append(0)
                mo_est_rng.append(0)
            mo_act.append(mo)

        # estimating mu peak
        if len(spikes) > 2:
            isi = diff(spikes)
            minisi = min(isi)
            minisi_index = argmin(isi)
            min_iii = range((spikes[minisi_index]/dt).astype(int),
                    (spikes[minisi_index+1]/dt).astype(int))
            min_iii = array(min_iii)+1
            K_min = minisi/dt
            mpe = vth/(tau*K_min*(1-exp(-dt/tau))) +\
                                            1/(tau*K_min)*sum(V[min_iii])*volt
            mp_est.append(mpe)

            # mu offs estimation as part of mu peak est
            if mo_est_rng[-1]:
                base_isi = get_misi_segments(d, 1, 4)
           #base_isi_index = argmin(abs(base_isi-isi))
           #base_iii = range((spikes[base_isi_index]/dt).astype(int),
           #        (spikes[base_isi_index+1]/dt).astype(int))
           #base_iii = array(base_iii)+1
           #K_base = base_isi/dt
           #moe = vth/(tau*K_base*(1-exp(-dt/tau))) +\
           #                                1/(tau*K_base)*sum(V[base_iii])*volt
                mo_est.append(vth/base_isi)
            else:
                mo_est.append(0)
        else:
            mp_est.append(0)
            mo_est.append(0)
        mp_act.append(mp)

    mp_act = array(mp_act)
    mp_est = array(mp_est)
    mo_est_rng = array(mo_est_rng)
    mo_est = array(mo_est)
    mo_act = array(mo_act)
    ma_act = array(ma_act)
    mw_ratios = array(mw_ratios)


'''
Reconstruct the input signal mu based on the Bibbona estimator
'''
def bibbona_reconstruction(d):
    tau = 10*ms
    vth = 10*mV
    dt = 0.1*ms
    spikes = d.get('spikes')[0]
    V = d.get('mem')[0]
    duration = d.get('duration')
    ISIs = diff(spikes)
    recon = zeros(duration/dt)
    for T, sp in zip(ISIs, spikes[:-1]):
        K = T/dt
        iii = (arange(sp*second, (sp+T)*second, dt)/dt+1).astype(int)
        m_est = vth/(tau*K*(1-exp(-dt/tau))) + sum(V[iii])*volt/(tau*K)
        recon[sp/dt] = m_est

    return recon


'''
Calculate the mean slope of the pre-spike membrane potential for all spikes
in each part of the sine wave period, after folding.
'''
def slope_bin_avg(slopes, spikes, freq, binwidth):
    '''
    Calculates the average value of the slopes within each bin of size
    binwidth after the slopes are aligned to the sine wave period.
    '''
    if not len(slopes):
        return zeros(1/(freq*binwidth))
    if is_dimensionless(freq):
        freq *= Hz
    if is_dimensionless(binwidth):
        binwidth *= second
    period = 1/freq
    periodspikes = spikes % period
    bins = frange(0*second, period, binwidth)
    avgslopes = []
    for l, r in zip(bins[:-1], bins[1:]):
        binspike_inds = bitwise_and(periodspikes > l, periodspikes < r)
        binslopes = slopes[binspike_inds]
        avgslopes.append(mean(binslopes) if len(binslopes) else 0)
    return array(avgslopes)

'''
Plot estimated vs actual with error bars
'''
def plot_est_vs_act_err(x_data, y_data, x_label, y_label, figname='',
                                                        plot_zeros=False):
    rc('text', usetex=True) # use latex
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
    xvalues([5, 10, 15, 20], fontsize=fontsize)
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
    #xvalues(fontsize=fontsize)
    #axis(xmin=0, xmax=3)
    #filename_pdf = "%s/%s_err.pdf" % (plot_save_loc, figname)
    #filename_png = "%s/%s_err.png" % (plot_save_loc, figname)
    #savefig(filename_pdf)
    #savefig(filename_png)
    #print "Figures saved: %s and %s" % (filename_pdf, filename_png)
    #clf()
    #print "\n\n"
    return

def discrete_sigma_recon(spikes, mu_t_est, V, dt=0.1*ms):
    '''
    Estimates the value of sigma for each interspike interval based on the
    mu estimator from Bibbona et al., 2008.
    '''
    isi = diff(spikes)
    isi = insert(isi, 0, spikes[0])
    recon = []
    for T, sp, mu in zip(isi, spikes, mu_t_est):
        K = float(T/dt)
        iii = (arange((sp-T)*second, (sp)*second, dt)/dt).astype(int)
        thesum = sum(((V[iii]-mu*tau+(mu*tau-V[iii-1])*exp(-dt/tau))**2)/(tau*(1-exp(-2*dt/tau))))
        s_est = sqrt(2/(K-1) * thesum)
        recon.append(s_est)
    return array(recon)

def discrete_sigma_recon_avg(spikes, mu_t_est, V, freq, nbins, dt=0.1*ms):
    '''
    Estimates the value of sigma for each interspike interval based on the
    sigma estimator from Bibbona et al., 2008, folds the estimated values
    based on period using the frequency `freq` and averages the estimated
    values in `nbins` bins.

    TODO: Try a moving average instead of binning.
    '''

    recon = discrete_sigma_recon(spikes, mu_t_est, V, dt)
    if is_dimensionless(freq):
        freq *= Hz
    period = 1/freq
    periodspikes = spikes % period
    bins = linspace(0*second, period, nbins)
    avgsigma = []
    for l, r in zip(bins[:-1], bins[1:]):
        binspike_inds = bitwise_and(periodspikes > l, periodspikes < r)
        binsigma = recon[binspike_inds]
        avgsigma.append(mean(binsigma) if len(binsigma) else 0)
    return array(avgsigma)

def plot_errors_vs_feat(feat, archive, bins=None):
    '''
    Feature should be a string of the feature to plot against.
    Archive should be a dictionary containing all the features
                                                    (parameters and data).
    '''
    feature = archive[feat]
    if bins is None:
        xvalues = unique(feature)
    else:
        xvalues = bins
    archive_keys = archive.keys()
    actual_keys = [k for k in archive_keys if k.endswith('actual')]
    estimation_errors = {}
    figure()
    clf()
    hold(True)
    nplots = len(actual_keys)
    for n, act_key in enumerate(actual_keys):
        base_key = act_key.replace('_actual', '')
        est_key = base_key+'_est'
        actuals = archive[act_key]  # maybe catch an exception here
        estimates = archive[est_key]
        errors = [abs(a-e)/a if a else e for a, e in zip(actuals, estimates)]
        errors = array(errors)
        estimation_errors[base_key] = errors
        #subplot(nplots, 1, n+1)
        if len(xvalues) == len(feature):
            # no grouping required
            plot(xvalues, errors, '-o', label=base_key+' error')
        else:
            # requires error bars
            digibins = digitize(feature, xvalues)-1
            udb = unique(digibins)
            means = [mean(errors[digibins == i]) for i in udb]
            stds  = [ std(errors[digibins == i]) for i in udb]
            x = unique(xvalues[digibins])
            #errorbar(xvalues[:-1], means, yerr=stds, fmt='-o',
            #            label=base_key+' error')
            plot(x, means, '-o', label=base_key+' error')
        #ylabel(base_key+' error')
        xlabel(feat)
    legend(loc='best')
    return feature, estimation_errors


from itertools import combinations, product
def stochastic_resonance(archive):
    '''
    Finds the noise amplitude which maximises the difference in output rates
    of two input mu configurations (mu_offs and mu amp).
    '''
    ioff()
    mu_offs_actual = archive['mu_offs_actual']
    mu_amp_actual  = archive['mu_amp_actual']
    mu_peaks_actual = archive['mu_peaks_actual']
    sigma_offs_actual = archive['sigma_offs_actual']
    sigma_amp_actual = archive['sigma_amp_actual']
    sigma_peaks_actual = archive['sigma_peaks_actual']
    freq_actual = archive['freq_actual']
    freq_est = archive['freq_est']
    outrates = archive['outrates']
    mo_unique = unique(mu_offs_actual)
    ma_unique = unique(mu_amp_actual)
    so_unique = unique(sigma_offs_actual)
    sa_unique = unique(sigma_amp_actual)
    # stick to one frequency
    alldiffs = []
    small_errs = []
    large_errs = []
    for f in unique(freq_actual):
        f_inds = freq_actual == f
        for ma in ma_unique:
            ma_inds = mu_amp_actual == ma
            mo_pairs = combinations(mo_unique, 2)
            for apair in mo_pairs:
                small = min(apair)
                large = max(apair)
                # find sigma which produces biggest outrate difference
                small_inds = bitwise_and(mu_offs_actual == small, ma_inds)
                large_inds = bitwise_and(mu_offs_actual == large, ma_inds)
                small_inds = bitwise_and(small_inds, f_inds)
                large_inds = bitwise_and(large_inds, f_inds)
                pair_diffs = []
                sigma_confs = product(sa_unique, so_unique)
                for sa, so in sigma_confs:
                    # get the firing rates for each unique combo
                    sigma_inds = bitwise_and(sigma_amp_actual == sa,
                                                sigma_offs_actual == so)
                    small_sim_ind = bitwise_and(sigma_inds, small_inds)
                    large_sim_ind = bitwise_and(sigma_inds, large_inds)
                    small_sim_ind = flatnonzero(small_sim_ind)
                    large_sim_ind = flatnonzero(large_sim_ind)
                    if len(small_sim_ind) == len(large_sim_ind) > 0:
                        print("Firing rates for ma %f, so %f, sa %f, f %f" % (
                            ma, so, sa, f))
                        print("mo %f: %s" % (small, outrates[small_sim_ind]))
                        print("mo %f: %s" % (large, outrates[large_sim_ind]))
                        rate_diff = max(outrates[large_sim_ind])-\
                                min(outrates[small_sim_ind])
                        pair_diffs.append((sa, so, rate_diff))
                        alldiffs.append(rate_diff)
                        sm_err = abs(
                                freq_actual[small_sim_ind]-freq_est[small_sim_ind]
                                    )/freq_actual[small_sim_ind]
                        lg_err = abs(
                                freq_actual[large_sim_ind]-freq_est[large_sim_ind]
                                    )/freq_actual[large_sim_ind]
                        small_errs.append(sm_err)
                        large_errs.append(lg_err)
                maxdiff = max(pair_diffs, key=lambda rd: rd[2])
                print("Outrate difference maximized at sa %f, so %f. Diff: %f" % (
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
    figure()
    clf()
    plot(alldiffs, small_errs, 'o', label='small')
    plot(alldiffs, large_errs, 'o', label='large')
    legend(loc='best')
    xlabel('Outrate difference')
    ylabel('Freq. est. error')
    savefig('figures/stoch_res_diff/stochres_vs_fr_est_err.png')



