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

duration = 2*second
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
            #'mem': mem_mon.values,
            'spikes': st_mon.spiketimes,
            'mslope': mslope,
            'slopes': slopes
            }

if __name__=='__main__':
    data_dir = 'fullmo_rand'
    num_mo = 100; num_ma = 5; num_so = 5; num_sa = 5; num_freq = 5;
    mo_min = 10; mo_max = 60
    mu_offs = random(num_mo)*(mo_max-mo_min)+mo_min
    ma_min = 10; ma_max = 60
    mu_amp = random(num_ma)*(ma_max-ma_min)+ma_min
    so_min = 0; so_max = 60
    sigma_offs = random(num_so)*(so_max-so_min)+so_min
    sa_min = 0; sa_max = 60
    sigma_amp = random(num_sa)*(sa_max-sa_min)+sa_min
    freq_min = 1; freq_max = 30
    freq = random(num_freq)*(freq_max-freq_min)+freq_min
    configs = zip(mu_offs, mu_amp, sigma_offs, sigma_amp, freq)
    numsims = num_mo*num_ma*num_so*num_sa*num_freq
    data = DataManager(data_dir)
    print("Simulations configured. Running ...")
    run_tasks(data, ousim, configs, gui=False, verbose=False, numitems=numsims)
    print("Simulations done!")

