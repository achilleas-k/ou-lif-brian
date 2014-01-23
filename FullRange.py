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
N = 1   # number of simulations per configuration
dt = defaultclock.dt
tau = 20*ms
v_th = 20*mV
t_refr = 2*ms
v_reset = 0*mV

def ousim(mu_offs, mu_amp, sigma_offs, sigma_amp, freq, report):
    clear(True)
    gc.collect()
    reinit_default_clock()
    mu_offs = mu_offs*mV
    mu_amp = mu_amp*mV
    sigma_amp = sigma_amp*mV
    sigma_offs = sigma_offs*mV if sigma_offs*mV > sigma_amp else sigma_amp
    freq = freq*Hz
    freq_ang = freq*2*pi    # angluar frequency for equation
    eqs=Equations('dI/dt = -I/dt + xi*tau**-0.5 : 1')
    eqs+=Equations('dV/dt = (mu-V)/tau + sigma*I/dt : volt')
    eqs+=Equations('mu = mu_amp*sin(t*freq_ang) + mu_offs : volt')
    eqs+=Equations('dsigma/dt = sigma_amp*cos(t*freq_ang)*freq_ang : volt')
    eqs.prepare()
    group = NeuronGroup(N, eqs, threshold=v_th, refractory=t_refr,
            reset=v_reset)
    group.sigma = sigma_offs
    #mu_mon = StateMonitor(group, 'mu', record=True)
    #sigma_mon = StateMonitor(group, 'sigma', record=True)
    mem_mon = StateMonitor(group, 'V', record=True)
    st_mon = SpikeMonitor(group)
    run(duration, report=None)
    mem_mon.insert_spikes(st_mon, value=v_th)
    mslope = 0
    slopes = array([])
    for i in range(N):
        cmslope, cslopes = nt.firing_slope(mem_mon[i], st_mon[i])
        mslope += cmslope
        slopes = append(slopes, cslopes)
    mslope /= N
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
    data_dir = "fullrange"
    mu_offs = linspace(0, 40, 6)
    mu_amp = linspace(0, 40, 6)
    sigma_offs = linspace(0, 40, 6)
    sigma_amp = linspace(0, 40, 6)
    freq = [1, 5, 10, 15]
    configs = itertools.product(mu_offs, mu_amp, sigma_offs, sigma_amp, freq)
    numsims=len(mu_offs)*len(mu_amp)*len(sigma_offs)*len(sigma_amp)*len(freq)
    data = DataManager(data_dir)
    print("Simulations configured. Running ...")
    run_tasks(data, ousim, configs, gui=False, verbose=False, numitems=numsims)
    print("Simulations done!")

