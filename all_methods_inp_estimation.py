from brian import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import sys
import gc
import itertools
from time import time
import neurotools as nt
from scipy.stats import halfnorm, kstest

forever = 1e10*second
duration = 3*second
N = 10
dt = defaultclock.dt
v_th = 10*mV
tau = 10*ms
t_refr = 0*ms
v_reset = 0*mV
V0 = 0*mV
min_spikes = 2000


'''
Monte carlo sampling parameters
'''
nconfigs = 200
mu_offs_min, mu_offs_max = 0, 2
mu_amp_min, mu_amp_max = 0, 2
sigma_offs_min, sigma_offs_max = 0, 1
sigma_amp_min, sigma_amp_max = 0, 1
freq_min, freq_max = 1, 30

def ousim(mu_offs, mu_amp, sigma_offs, sigma_amp, freq, report):
    clear(True)
    gc.collect()
    reinit_default_clock()
    np.random.seed(int(time()+(mu_offs+mu_amp+sigma_offs+sigma_amp)*1e10))
    mu_offs = mu_offs*mV/ms
    mu_amp = mu_amp*mV/ms
    sigma_amp = sigma_amp*mV/sqrt(ms)
    sigma_offs = sigma_offs*mV/sqrt(ms)\
            if sigma_offs*mV/sqrt(ms) > sigma_amp else sigma_amp
    freq = freq*Hz
    freq_ang = freq*2*pi # angluar frequency for equation

    eqs =Equations('dV/dt = mu-(V+V0)/tau + sigma*I/sqrt(dt) : volt')
    eqs+=Equations('dI/dt = -I/dt + xi/sqrt(dt) : 1')
    eqs+=Equations('mu = mu_amp*sin(t*freq_ang) + mu_offs : volt/second')
    eqs+=Equations('sigma = sigma_amp*sin(t*freq_ang) + sigma_offs :'
                                                        ' volt/sqrt(second)')
    eqs.prepare()
    group = NeuronGroup(N, eqs, threshold=v_th, refractory=t_refr,
            reset=v_reset)
    group.V = V0
    mu_mon = StateMonitor(group, 'mu', record=True)
    sigma_mon = StateMonitor(group, 'sigma', record=True)
    I_mon = StateMonitor(group, 'I', record=True)
    mem_mon = StateMonitor(group, 'V', record=True)
    st_mon = SpikeMonitor(group)

    #if (mu_amp+mu_offs)*tau < v_th:
    #    # subthreshold: might take too long to reach target number of
    #    # spikes; reduce duration
    #    duration = 2*second
    #else:
    #    duration = forever
    #duration = forever
    #@network_operation
    #def stop_condition():
    #    if defaultclock.t > 1*second and st_mon.nspikes < 5:
    #        stop()
    #    if st_mon.nspikes >= min_spikes:
    #        stop()

    run(duration, report=report)
    mem_mon.insert_spikes(st_mon, value=v_th)
    mslope = 0
    slopes = []
    for i in range(N):
        cmslope, cslopes = nt.firing_slope(mem_mon[i], st_mon[i])
        mslope += cmslope
        slopes.extend(cslopes)
    mslope /= N
    print("Sim with MO %s MA %s and SO %s "
            "SA %s and freq %s fired %i spikes in %s\n\n" % (
            display_in_unit(mu_offs, mV/ms),
            display_in_unit(mu_amp, mV/ms),
            display_in_unit(sigma_offs, mV/sqrt(ms)),
            display_in_unit(sigma_amp, mV/sqrt(ms)),
            display_in_unit(freq, Hz),
            st_mon.nspikes,
            defaultclock.t))

    return {
            'mu_offs': mu_offs,
            'mu_amp': mu_amp,
            'sigma_offs': sigma_offs,
            'sigma_amp': sigma_amp,
            'freq': freq,
            #'mem': mem_mon.values,
            'spikes': st_mon.spiketimes,
            'mslope': mslope,
            'slopes': slopes,
            #'mu': mu_mon.values,
            #'sigma': sigma_mon.values,
            #'xi': I_mon.values,
            }

def sinewave(offs, amp, freq, t):
    return offs+amp*sin(2*pi*freq*t)

def plot_est_vs_act(est, actual, x_label, y_label):
    good_thr = 0.1
    est = array(est)
    actual = array(actual)
    relerr = abs(est-actual)/actual
    good_ratio = 1.0*count_nonzero(relerr < good_thr)/len(relerr)
    oneline = linspace(0, max(actual)+max(actual)/10)
    figure(num=1, figsize=(13.66,7.68), dpi=100)
    plot(est, actual, marker='.', linestyle='', markersize=11,
            color='blue', label=x_label)
    plot(oneline, oneline, color='red', linestyle='--', label='1:1')
    xlabel(x_label)
    ylabel(y_label)
    legend(loc='best')
    title('%s vs %s. %f %% have an error < %0.2f '
            'Mean error: %f' % (
        y_label, x_label,
        good_ratio*100,
        good_thr,
        mean(relerr)))
    figname = "tmp/%s_vs_%s.png" % (y_label, x_label)
    savefig(figname)
    print "Figure saved: %s" % figname
    clf()
    return

def fft_binst_est(st, dt):
    '''
    Estimates for the frequency and mu parameters based on the fourier
    transform of the binned spike train
    '''
    binst = nt.times_to_bin(st)
    fft_st = fft(binst)
    cutfreq = len(fft_st)/2
    fft_st[cutfreq:] = zeros(len(fft_st[cutfreq:]))
    base_amp = fft_st[0]
    fft_st[0] = 0
    max_amp_i = int(where(fft_st == max(fft_st))[0])
    max_freq = max_amp_i/(len(binst)*dt)
    max_amp = fft_st[max_amp_i]
    return base_amp, max_amp, max_freq

def fft_acorr_est(st, dt):
    '''
    Estimates for the frequency and mu parameters based on the fourier
    transform of the autocorrelogram of the spike train
    '''
    width = st[-1]/2
    bin = 10*ms
    acorr = autocorrelogram(st, width=width, bin=bin)
    acorr = acorr[len(acorr)/2:] # acorr is symmetric
    fft_acorr = fft(acorr)
    cutfreq = len(fft_acorr)/2
    fft_acorr[cutfreq:] = zeros(len(fft_acorr[cutfreq:])) # low pass filter
    base_amp = fft_acorr[0]*dt
    fft_acorr[0] = 0
    max_amp_i = int(where(fft_acorr == max(fft_acorr))[0])
    max_freq = max_amp_i/(len(fft_acorr)*bin)
    max_amp = fft_acorr[max_amp_i]*dt
    est_freq = max_freq
    return base_amp, max_amp, max_freq

def mean_rate_at_peak(spikes, freq):
    num_periods = round(sum([max(st) for st in spikes.itervalues()\
            if len(st) > 0])*second*freq)
    if not num_periods:
        return 0
    period_spikes = []
    period = 1./freq
    for st in spikes.itervalues():
        period_spikes.extend(mod(st, period))
    period_peak_start = float(period*1./8)
    period_peak_end = float(period*3./8)
    peak_spikes = [t for t in period_spikes if t > period_peak_start and\
            t < period_peak_end]
    spikes_per_peak = len(peak_spikes)/num_periods
    rate_at_peak = spikes_per_peak/(period_peak_end-period_peak_start)
    return rate_at_peak

if __name__ == '__main__':
    data_dir = sys.argv[1]
    data = DataManager(data_dir)
    print("\n")
    if '--no-run' not in sys.argv:
        mu_offs = mu_offs_min+rand(nconfigs)*mu_offs_max
        mu_amp = mu_amp_min+rand(nconfigs)*mu_amp_max
        sigma_offs = sigma_offs_min+rand(nconfigs)*sigma_offs_max
        sigma_amp = sigma_amp_min+rand(nconfigs)*sigma_amp_max
        freq = freq_min+rand(nconfigs)*freq_max
        #nsims=len(mu_amp)*len(mu_offs)*len(sigma_amp)*len(sigma_offs)*len(freq)
        #params_prod = itertools.product(mu_offs, mu_amp, sigma_offs, sigma_amp,
        #        freq)
        params = zip(mu_offs, mu_amp, sigma_offs, sigma_amp, freq)
        nsims = nconfigs
        print "Created configuration generator"
        print "Simulations configured. Running ..."
        run_tasks(data, ousim, params, gui=True, poolsize=4,
                numitems=nsims)
        print "Simulations done!"
    else:
        print "Skipping simulation run. Working with %s.data directory\n" %\
                data_dir

    numsims = data.itemcount()
    print "Total number of simulations: %i" % numsims

    mu_offs_est_a = []
    mu_offs_est_am = []
    mu_offs_est_b = []
    mu_offs_est_bm = []
    mu_offs_actual = []
    mu_offs_actual10 = []
    mu_amp_est_a = []
    mu_amp_est_am = []
    mu_amp_est_b = []
    mu_amp_est_bm = []
    mu_amp_actual = []
    mu_amp_actual10 = []
    mu_peaks = []
    mu_means = []
    mu_est = []
    mu_peak_est_mrap = []
    sigma_est = []
    sigma_means = []
    sigma_peaks = []
    freq_actual = []
    freq_actual10 = []
    freq_est_a = []
    freq_est_am = []
    freq_est_b = []
    freq_est_bm = []
    frate_peak = []

    i = 1
    for d in data.itervalues():
        print "\rProcessing %i of %i" % (i, numsims),
        sys.stdout.flush()
        slopes = divide(d.get('slopes'), 1000)
        if len(slopes) == 0:
            print "No data available (o) --  skipping."
            i+=1
            continue

        spikes = d.get('spikes')
        mu_offs = d.get('mu_offs')
        mu_amp = d.get('mu_amp')
        sigma_offs = d.get('sigma_offs')
        sigma_amp = d.get('sigma_amp')
        freq = d.get('freq')
        filehead = "tmp/m%f-%fs%f-%f" % (mu_offs/(mV/ms),
                mu_amp/(mV/ms),
                sigma_offs/(mV/sqrt(ms)),
                sigma_amp/(mV/sqrt(ms)))

        mu_p = mu_amp+mu_offs
        sigma_p = sigma_amp+sigma_offs

        firing_mu = []
        firing_sigma = []
        base_amps_a = []
        max_amps_a = []
        max_freqs_a = []
        base_amps_b = []
        max_amps_b = []
        max_freqs_b = []
        isis = []
        for st in spikes.values():
            if len(st) < 10:
                continue

            firing_mu.extend(sinewave(mu_offs, mu_amp, freq, st))
            firing_sigma.extend(sinewave(sigma_offs, sigma_amp, freq, st))

            base_amp_b, max_amp_b, max_freq_b = fft_binst_est(st, dt)
            base_amp_a, max_amp_a, max_freq_a = fft_acorr_est(st, dt)
            base_amps_a.append(base_amp_a)
            max_amps_a.append(max_amp_a)
            max_freqs_a.append(max_freq_a)
            base_amps_b.append(base_amp_b)
            max_amps_b.append(max_amp_b)
            max_freqs_b.append(max_freq_b)
            isis.extend(diff(st))

        if not firing_mu and not firing_sigma:
            'No data in any of the spike trains'
            i+=1
            continue

        h, b = histogram(isis, bins=100, density=True)
        frate_peak.append(b[h == max(h)][0])
        mu_mean_est = v_th / array(frate_peak)

        mean_base_amp_a = mean(base_amps_a)
        mean_max_amp_a = mean(max_amps_a)
        mean_freq_a = mean(max_freqs_a)
        mean_base_amp_b = mean(base_amps_b)
        mean_max_amp_b = mean(max_amps_b)
        mean_freq_b = mean(max_freqs_b)


        mu_offs_actual.append(mu_offs)
        mu_amp_actual.append(mu_amp)
        freq_actual.append(freq)

        mu_offs_actual10.extend(ones(len(base_amps_a))*mu_offs)
        mu_amp_actual10.extend(ones(len(max_amps_a))*mu_amp)
        freq_actual10.extend(ones(len(max_freqs_a))*freq)

        mu_offs_est_a.extend(base_amps_a) # mu_offs estimates (acorr)
        mu_amp_est_a.extend(max_amps_a)  # mu_amp estimates (acorr)
        freq_est_a.extend(max_freqs_a) # frequency estimates (acorr)

        mu_offs_est_am.append(mean_base_amp_a)
        mu_amp_est_am.append(mean_max_amp_a)
        freq_est_am.append(mean_freq_a)

        mu_offs_est_b.extend(base_amps_b) # mu_offs estimates (bin)
        mu_amp_est_b.extend(max_amps_b)  # mu_amp estimates (bin)
        freq_est_b.extend(max_freqs_b) # frequency estimates (bin)

        mu_offs_est_bm.append(mean_base_amp_b)
        mu_amp_est_bm.append(mean_max_amp_b)
        freq_est_bm.append(mean_freq_b)

        mrap = mean_rate_at_peak(spikes, mean_freq_a)

        mu_peak_est_mrap.append(multiply(v_th, mrap)+0.6)
        mu_m = mean(firing_mu)
        sigma_m = mean(firing_sigma)

        mu_peaks.append(mu_p)
        mu_means.append(mu_m)
        sigma_peaks.append(sigma_p)
        sigma_means.append(sigma_m)

        mslope = mean(slopes)
        vslope = var(slopes)

        Vn_mean = mslope**2 * pi/2 # estimation of Vn based on mean
        sigma_est.append(sqrt(2*Vn_mean / (tau * (1+exp(-2*dt/tau)))))
                                                        # sigma estimates
        i += 1

    plot_est_vs_act(freq_est_a, freq_actual10, 'freq_est acorr', 'freq_act')
    plot_est_vs_act(freq_est_b, freq_actual10, 'freq_est bin', 'freq_act')
    plot_est_vs_act(freq_est_am, freq_actual,
                                    'freq_est acorr (mean)', 'freq_act')
    plot_est_vs_act(freq_est_bm, freq_actual,
                                    'freq_est bin (mean)', 'freq_act')

    plot_est_vs_act(sigma_est, sigma_peaks, 'sigma_est', 'sigma_peaks')
    plot_est_vs_act(sigma_est, sigma_means, 'sigma_est', 'sigma_means')

    plot_est_vs_act(mu_mean_est, mu_means, 'mu_mean_est', 'mu_means')
    plot_est_vs_act(mu_offs_est_am, mu_offs_actual,
                            'mu_offs est acorr (mean)', 'mu_offs act')
    plot_est_vs_act(divide(mu_offs_est_bm, 200), mu_offs_actual,
                                'mu_offs est bin (mean)', 'mu_offs act')

    plot_est_vs_act(mu_offs_est_a, mu_offs_actual10,
                                    'mu_offs est acorr', 'mu_offs act')
    plot_est_vs_act(divide(mu_offs_est_b, 200), mu_offs_actual10,
                                        'mu_offs est bin', 'mu_offs act')

    plot_est_vs_act(mu_peak_est_mrap, mu_peaks, 'mu_peak est mrap', 'mu_peaks')

    from IPython import embed
    embed()

