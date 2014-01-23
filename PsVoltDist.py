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
N = 10
dt = defaultclock.dt
v_th = 10*mV
tau = 10*ms
t_refr = 0*ms
v_reset = 0*mV
V0 = 0*mV
min_spikes = 2000


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
    eqs+=Equations('sigma = sigma_amp*sin(t*freq_ang) + sigma_offs : volt/sqrt(second)')
    eqs.prepare()
    group = NeuronGroup(N, eqs, threshold=v_th, refractory=t_refr,
            reset=v_reset)
    group.V = V0
    mu_mon = StateMonitor(group, 'mu', record=True)
    sigma_mon = StateMonitor(group, 'sigma', record=True)
    I_mon = StateMonitor(group, 'I', record=True)
    mem_mon = StateMonitor(group, 'V', record=True)
    st_mon = SpikeMonitor(group)

    if (mu_amp+mu_offs)*tau < v_th:
        # subthreshold: might take too long to reach 1k spikes
        # reduce duration
        duration = 5*second
    else:
        duration = forever

    @network_operation
    def stop_condition():
        if st_mon.nspikes >= min_spikes:
            stop()

    run(duration, report=report)
    mem_mon.insert_spikes(st_mon, value=v_th)
    mslope = 0
    slopes = array([])
    for i in range(N):
        cmslope, cslopes = nt.firing_slope(mem_mon[i], st_mon[i])
        mslope += cmslope
        slopes = append(slopes, cslopes)
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
            'mem': mem_mon.values,
            'spikes': st_mon.spiketimes,
            'mslope': mslope,
            'slopes': slopes,
            #'mu': mu_mon.values,
            #'sigma': sigma_mon.values,
            #'xi': I_mon.values,
            }


if __name__ == '__main__':
    data_dir = sys.argv[1]
    data = DataManager(data_dir)
    print("\n\n")
    if '--no-run' not in sys.argv:
        mu_offs = [1]
        mu_amp = [0.5]
        sigma_offs = [0.3]
        sigma_amp = [0.2]
        freq = [15]
        nsims=len(mu_amp)*len(mu_offs)*len(sigma_amp)*len(sigma_offs)*len(freq)
        params_prod = itertools.product(mu_offs, mu_amp, sigma_offs, sigma_amp,
                freq)
        print "Created configuration generator"
        print("Simulations configured. Running ...")
        run_tasks(data, ousim, params_prod, gui=True, poolsize=4,
                numitems=nsims)
        print("Simulations done!")
    else:
        print "Skipping simulation run. Working with %s.data directory" %\
                data_dir

    numsims = data.itemcount()
    print("Total number of simulations: %i" % numsims)

    rms_mean_results = array([])
    rms_peak_results = array([])

    i = 1
    for d in data.itervalues():
        print "\rProcessing %i of %i" % (i, numsims),
        sys.stdout.flush()
        slopes = d.get('slopes')
        if len(slopes) == 0:
            print "(o): No data available -- skipping %i." % i
            i+=1
            continue

        spikes = d.get('spikes')
        mu_offs = d.get('mu_offs')
        mu_amp = d.get('mu_amp')
        sigma_offs = d.get('sigma_offs')
        sigma_amp = d.get('sigma_amp')
        freq = d.get('freq')
        filehead = 'm%f-%fs%f-%f' % (mu_offs/(mV/ms),
                mu_amp/(mV/ms),
                sigma_offs/(mV/sqrt(ms)),
                sigma_amp/(mV/sqrt(ms)))

        def sinewave(offs, amp, freq, t):
            return offs+amp*sin(2*pi*freq*t)

        firing_mu = []
        firing_sigma = []
        for st in spikes.values():
            firing_mu.extend(sinewave(mu_offs, mu_amp, freq, st))
            firing_sigma.extend(sinewave(sigma_offs, sigma_amp, freq, st))

        psvolt = v_th - array(slopes)*dt
        pdf_actual, edges = histogram(psvolt, density=True, bins=20)
        dx = diff(edges)[0]
        x = array(edges[:-1]+dx/2) # x = bin centres

        mean_actual = mean(psvolt)
        var_actual = var(psvolt)

        mu_peak = mu_amp+mu_offs
        sigma_peak = sigma_amp+sigma_offs
        mu_mean = mean(firing_mu)*volt/second
        sigma_mean = mean(firing_sigma)*volt/sqrt(second)
        t = dt

        #Mt_peak = mu_peak*tau + (v_th - mu_peak*tau)*exp(-t/tau)
        Mt_peak = v_th
        Vt_peak = (((sigma_peak**2) * tau)/2) * (1-exp((-2*t)/tau))
        if Vt_peak == 0:
            i+=1
            continue
        pdf_peak = ((2*pi*Vt_peak)**-0.5) * exp(-((x-Mt_peak)**2) / (2*Vt_peak))
        pdf_peak /= trapz(pdf_peak, dx=dx) # rescale to 1
        mean_peak = average(x, weights=pdf_peak)
        var_peak = dot(pdf_peak, (x-mean_peak)**2)/pdf_peak.sum()

        #Mt_mean = mu_mean*tau + (v_th - mu_mean*tau)*exp(-t/tau)
        Mt_mean = v_th
        Vt_mean = (((sigma_mean**2) * tau)/2) * (1-exp((-2*t)/tau))
        pdf_mean = ((2*pi*Vt_mean)**-0.5) * exp(-((x-Mt_mean)**2) / (2*Vt_mean))
        pdf_mean /= trapz(pdf_mean, dx=dx) # rescale to 1
        mean_mean = average(x, weights=pdf_mean)
        var_mean = dot(pdf_mean, (x-mean_mean)**2)/pdf_mean.sum()

        f1 = figure(num=1, figsize=(13.66,7.68), dpi=100)
        plot(x, pdf_actual, color='blue', label='Numerical')
        plot(x, pdf_peak, color='red', linestyle='--', label='Theoretical')
        plot(x, pdf_mean, color='green', linestyle='--', label='Theo mean')
        xlabel('Membrane potential (volt)')
        title('Theoretical and numerical pre-spike volt distributions')
        legend(loc='best')
        figname = 'tmp/%s_psvolt.png' % filehead
        savefig(figname)
        f1.clf()
        rms_mean = sqrt(mean((pdf_mean - pdf_actual)**2))
        rms_peak = sqrt(mean((pdf_peak - pdf_actual)**2))

        hn_mean = halfnorm(loc=v_th-Mt_mean, scale=sqrt(Vt_mean)/dt)
        hn_peak = halfnorm(loc=v_th-Mt_peak, scale=sqrt(Vt_peak)/dt)

        quantiles = frange(0, 1, 0.01)

        hn_mean_qs = hn_mean.ppf(quantiles)
        hn_peak_qs = hn_peak.ppf(quantiles)

        #from IPython import embed
        #embed()
        #sys.exit()

        #print "RMS mean dist: %s\nRMS peak dist: %f" % (rms_mean, rms_peak)
        #print "KS Mean: %f\nKS Peak: %f" % (ks_mean[1], ks_peak[1])
        #print "Avgs Actual, Mean, Peak: %10.9f, %10.9f, %10.9f" % (
        #        mean_actual,
        #        mean_mean,
        #        mean_peak)
        #print "Vars Actual, Mean, Peak: %10.9f, %10.9f, %10.9f" % (
        #        var_actual,
        #        var_mean,
        #        var_peak)
        #print "ThAv Actual, Mean, Peak: %10.9f, %10.9f, %10.9f" % (
        #        0,
        #        v_th-sqrt(2*Vt_mean/pi),
        #        v_th-sqrt(2*Vt_peak/pi))
        #print "ThVr Actual, Mean, Peak: %10.9f, %10.9f, %10.9f" % (
        #        0,
        #        Vt_mean*(1-2/pi),
        #        Vt_peak*(1-2/pi))

        rms_mean_results = append(rms_mean_results, rms_mean)
        rms_peak_results = append(rms_peak_results, rms_peak)


        percentiles = list(quantiles*100) # percentile function fails with array
        actual_qs = percentile(slopes, q=percentiles)

        oneline = [0, max(actual_qs)]

        figure(num=2, figsize=(13.66,7.68), dpi=100)
        plot(hn_mean_qs, actual_qs, marker='.', linestyle='', markersize=15,
                color='green', label='Mean-based quantiles')
        plot(hn_peak_qs, actual_qs, marker='.', linestyle='', markersize=13,
                color='blue')
        plot(oneline, oneline, color='red', linestyle='--', label='1:1')
        xlabel('Theoretical quantiles')
        ylabel('Numerical quantiles')
        legend(loc='best')
        title('Q-Q plot')

        figname = 'tmp/%s_quant.png' % filehead
        savefig(figname)
        clf()
        i += 1

    figure(num=3, figsize=(13.66,7.68), dpi=100)
    plot(rms_mean_results, color='blue', marker='.',
            linestyle='', label='RMS (mean)')
    plot(movavg(rms_mean_results, 30), color='cyan', marker='',
            linestyle='-', linewidth=3, label='mov avg RMS (mean)')
    plot(rms_peak_results, color='red', marker='.',
            linestyle='', label='RMS (peak)')
    plot(movavg(rms_peak_results, 30), color='magenta', marker='',
            linestyle='-', linewidth=3, label='mov avg RMS (peak)')
    for i in range(len(rms_mean_results)):
        plot([i, i], [rms_mean_results[i], rms_peak_results[i]],
                color='blue' if rms_mean_results[i] < rms_peak_results[i]
                    else 'red')
    #yscale('log')
    xlabel('Sim # (arbitrary)')
    ylabel('RMS')
    legend(loc='best')
    title("RMS between actual and both mean and peak estimates."
            "\nMean wins: %i, Peak wins: %i" % (
                count_nonzero(rms_mean_results < rms_peak_results),
                count_nonzero(rms_mean_results > rms_peak_results)))
    savefig('tmp/rmsvalues.png')
    clf()
    #from IPython import embed
    #embed()
