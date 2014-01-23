import brian_no_units
from brian import *
from brian.library.random_processes import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import itertools
import neurotools as nt
import sys
import os
from datetime import datetime
import gc

duration = 1*second
N = 1 # number of simulations per configuration
dt = defaultclock.dt
tau = 20*ms
v_th = 20*mV
t_refr = 2*ms
v_reset = 0*mV
V0 = 0*mV

def ousim(mu_offs, mu_amp, sigma_offs, sigma_amp, freq, report):
    clear(True)
    gc.collect()
    reinit_default_clock()
    mu_offs = mu_offs*mV
    mu_amp = mu_amp*mV
    sigma_amp = sigma_amp*mV
    freq = freq*Hz
    freq_ang = freq*2*pi # angluar frequency for equation
    sigma_offs = sigma_offs*mV if sigma_offs*mV > sigma_amp else sigma_amp
    eqs=Equations('dV/dt = (mu-V+V0)/tau + sigma*I/dt : volt')
    eqs+=Equations('mu = mu_amp*sin(t*freq_ang) + mu_offs : volt')
    eqs+=Equations('dI/dt = -I/dt + xi*tau**-0.5 : 1')
    eqs+=Equations('dsigma/dt = sigma_amp*cos(t*freq_ang)*freq_ang : volt')
    eqs.prepare()
    group = NeuronGroup(N, eqs, threshold=v_th, refractory=t_refr,
            reset=v_reset)
    group.sigma = sigma_offs
    group.V = V0
    #mu_mon = StateMonitor(group, 'mu', record=True)
    #sigma_mon = StateMonitor(group, 'sigma', record=True)
    mem_mon = StateMonitor(group, 'V', record=True)
    st_mon = SpikeMonitor(group)
    run(duration, report=report)
    mem_mon.insert_spikes(st_mon, value=v_th)
    mslope = 0
    slopes = array([])
    for i in range(N):
        cmslope, cslopes = nt.firing_slope(mem_mon[i], st_mon[i])
        mslope += cmslope
        slopes = append(slopes, cslopes)
    mslope /= N
    print "Max slope for mo %f, ma %f, so %f, sa %f: %f" % (mu_offs, mu_amp, sigma_offs, sigma_amp, max(slopes))
    return {
            'mu_offs': mu_offs,
            'mu_amp': mu_amp,
            'sigma_offs': sigma_offs,
            'sigma_amp': sigma_amp,
            'freq': freq,
            'mem': mem_mon.values,
            'spikes': st_mon.spiketimes,
            'mslope': mslope,
            'slopes': slopes
            }


if __name__=='__main__':
    data_dir = 'findmax'
    mu_amp = linspace(1000, 2000, 2)
    mu_offs = linspace(1000, 2000, 2)
    sigma_amp = linspace(1000, 2000, 2)
    sigma_offs = linspace(1000, 2000, 2)
    freq = linspace(1, 100, 4)
    print "Created individual iterators"
    params_prod = itertools.product(mu_offs, mu_amp, sigma_offs, sigma_amp,
            freq)
    print "Created configuration generator (product)"
    data = DataManager(data_dir)
    print("Simulations configured. Running ...")
    run_tasks(data, ousim, params_prod, gui=False)
    print("Simulations done!")
