#import brian_no_units
from brian import *
from brian.library.random_processes import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import itertools
import neurotools as nt
import sys
import os
from time import time
from datetime import datetime
import gc

duration = 2*second
N = 1 # number of simulations per configuration-process
dt = defaultclock.dt
tau = 20*ms
v_th = 15*mV
t_refr = 0*ms
v_reset = 0*mV
V0 = 0*mV


def ousim(mu_offs, mu_amp, sigma_offs, sigma_amp, freq, report):
    clear(True)
    gc.collect()
    reinit_default_clock()
    np.random.seed(int(time()+(mu_offs+mu_amp+sigma_offs+sigma_amp)*1e10))
    mu_offs = mu_offs*mvolt/ms
    mu_amp = mu_amp*mvolt/ms
    sigma_amp = sigma_amp*mvolt/ms**0.5
    sigma_offs = sigma_offs*mvolt/ms**0.5\
            if sigma_offs*mvolt/ms**0.5 > sigma_amp else sigma_amp
    freq = freq*Hz
    freq_ang = freq*2*pi # angluar frequency for equation

    mu = mu_offs
    sigma = sigma_offs
    eqs=Equations('dV/dt = mu-(V+V0)/tau + sigma*xi : volt')
    #eqs+=Equations('dI/dt = -I/dt + xi*dt**-0.5 : 1')
    #eqs+=Equations('mu = mu_amp*sin(t*freq_ang)+mu_offs : volt/second')
    #eqs+=Equations('sigma = sigma_amp*sin(t*freq_ang)+sigma_offs : volt/second')
    eqs.prepare()
    group = NeuronGroup(N, eqs, threshold=v_th, refractory=t_refr,
            reset=v_reset)
    group.V = V0
    #mu_mon = StateMonitor(group, 'mu', record=True)
    #sigma_mon = StateMonitor(group, 'sigma', record=True)
    mem_mon = StateMonitor(group, 'V', record=True)
    st_mon = SpikeMonitor(group)
    run(duration, report=report)
    #subplot(311)
    #plot(mu_mon.times, mu_mon[0])
    #subplot(312)
    #plot(sigma_mon.times, sigma_mon[0])
    plot(mem_mon.times, mem_mon[0])
    suptitle("mu: %s, sigma: %s, nspikes: %i" % (mu_offs, sigma_offs, st_mon.nspikes))
    savefig("tmp/m%fs%f.png" % (mu_offs, sigma_offs))
    clf()
    #mem_mon.insert_spikes(st_mon, value=v_th)
    #mslope = 0
    #slopes = array([])
    #for i in range(N):
    #    cmslope, cslopes = nt.firing_slope(mem_mon[i], st_mon[i], w=2*ms)
    #    mslope += cmslope
    #    slopes = append(slopes, cslopes)
    #mslope /= N
    print "Sim with mu %s and sigma %s fired %i spikes" % (mu_offs,
            sigma_offs, st_mon.nspikes)
    return {
            'mu_offs': mu_offs,
            'mu_amp': mu_amp,
            'sigma_offs': sigma_offs,
            'sigma_amp': sigma_amp,
            'freq': freq,
            'mem': mem_mon.values,
            'spikes': st_mon.spiketimes,
            #'mslope': mslope,
            #'slopes': slopes
            }


if __name__=='__main__':
    data_dir = 'longdura_offs'
    mu_amp = linspace(0, 0, 1)
    mu_offs = linspace(0, 1, 10)[1:]
    sigma_amp = linspace(0, 0, 1)
    sigma_offs = linspace(0, 1, 10)[1:]
    freq = [10]
    nsims = len(mu_amp)*len(mu_offs)*len(sigma_amp)*len(sigma_offs)*len(freq)
    print "Created individual iterators"
    params_prod = itertools.product(mu_offs, mu_amp, sigma_offs, sigma_amp,
            freq)
    print "Created configuration generator (product)"
    data = DataManager(data_dir)
    print("Simulations configured. Running ...")
    run_tasks(data, ousim, params_prod, gui=True, poolsize=0)
    print("Simulations done!")
