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
min_spikes = 1000


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
        mu_offs = [0.2, 0.3, 0.5]
        mu_amp = [0.2, 0.3, 0.5]
        sigma_offs = [2]
        sigma_amp = [2]
        freq = [10]
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
    mu_est = array([])
    mu_actual = array([])
    sigma_est = array([])
    sigma_actual = array([])
    theta_est = array([])
    theta_actual = array([])
    i = 1
    for d in data.itervalues():
        print "\nProcessing %i of %i" % (i, numsims)
        slopes = d.get('slopes')
        if len(slopes) == 0:
            print "No data available (o) --  skipping."
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

        def sinewave(offs, amp, freq, t):
            return offs+amp*sin(2*pi*freq*t)

        mu_peak = mu_amp+mu_offs
        sigma_peak = sigma_amp+sigma_offs

        firing_mu = []
        firing_sigma = []
        for st in spikes.values():
            firing_mu.extend(sinewave(mu_offs, mu_amp, freq, st))
            firing_sigma.extend(sinewave(sigma_offs, sigma_amp, freq, st))


        slope_pdf, edges = histogram(slopes, density=True, bins=40)
        dx = diff(edges)[0]
        x = array(edges[:-1]+dx/2) # x = bin centres

        mslope = mean(slopes)
        vslope = var(slopes)

        Vn = vslope/(1-2/pi)
        sigma_est = append(sigma_est, sqrt(2*Vn / (tau * (1+exp(-2*dt/tau)))))

        sigma_actual = append(sigma_actual, sigma_peak)


        i += 1

    figure(figsize=(13.66,7.68), dpi=100)
    plot(sigma_est, sigma_actual, marker='.', linestyle='', markersize=15)
    savefig("tmp/sigma_estimation.png")


