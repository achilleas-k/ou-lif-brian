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

duration = 1.0*second
N = 1 # number of simulations per configuration
dt = defaultclock.dt
tau = 10*ms
v_th = 15*mV
t_refr = 0*ms
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
    eqs+=Equations('dI/dt = -I/dt + xi*tau**-0.5 : 1')
    eqs+=Equations('mu = mu_amp*sin(t*freq_ang) + mu_offs : volt')
    eqs+=Equations('sigma = sigma_amp*sin(t*freq_ang) + sigma_offs : volt')
    eqs.prepare()
    group = NeuronGroup(N, eqs, threshold=v_th, refractory=t_refr,
            reset=v_reset)
    group.V = V0
    mu_mon = StateMonitor(group, 'mu', record=True)
    sigma_mon = StateMonitor(group, 'sigma', record=True)
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
            'mu': mu_mon.values,
            'sigma': sigma_mon.values,
            }

if __name__=='__main__':
    data_dir = 'changeme'
    ma_start = 0; ma_end = 61; ma_step = 5
    mu_amp = xrange(ma_start, ma_end, ma_step)
    mo_start = 0; mo_end = 61; mo_step = 5
    mu_offs = xrange(mo_start, mo_end, mo_step)
    sa_start = 0; sa_end = 1; sa_step = 5
    sigma_amp = xrange(sa_start, sa_end, sa_step)
    so_start = 0; so_end = 1; so_step = 5
    sigma_offs = xrange(so_start, so_end, so_step)
    freq = [1, 5, 10, 15]
    print "Created individual iterators"
    params_prod = itertools.product(mu_offs, mu_amp, sigma_offs, sigma_amp,
            freq)
    print "Created configuration generator (product)"
    data = DataManager(data_dir)
    numsims = ((ma_end-ma_start)/ma_step+1) *\
        ((sa_end-sa_start)/sa_step+1) *\
        len(freq) *\
        ((mo_end - mo_start)/mo_step+1) *\
        ((so_end - so_start)/so_step+1)
    print("Simulations configured. Running ...")
    run_tasks(data, ousim, params_prod, gui=False, numitems=numsims)
    print("Simulations done!")
