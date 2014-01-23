from brian import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import sys
import gc
from time import time
import itertools

N = 1
nsims = 500
duration = 1*second
defaultclock.dt = dt = 0.1*ms
V_th = 10*mV
tau = 10*ms
t_refr = 0*ms
v_reset = 0*mV
V0 = 0*mV

'Random parameter ranges'
mo_min = 0; mo_max = 2
ma_min = 0; ma_max = 2

so_min = 0; so_max = 0.4
sa_min = 0; sa_max = 0.4

freq_min = 5; freq_max = 20

def ousim(mu_offs, mu_amp, sigma_offs, sigma_amp, freq, report):
    clear(True)
    gc.collect()
    reinit_default_clock()
    seed = int(time()+(mu_offs+mu_amp+sigma_offs+sigma_amp)*1e3)
    np.random.seed(seed)
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
    group = NeuronGroup(N, eqs, threshold='V>V_th', refractory=t_refr,
            reset=v_reset)
    group.V = V0
    #mu_mon = StateMonitor(group, 'mu', record=True)
    #sigma_mon = StateMonitor(group, 'sigma', record=True)
    #I_mon = StateMonitor(group, 'I', record=True)
    mem_mon = StateMonitor(group, 'V', record=True)
    st_mon = SpikeMonitor(group)

    run(duration, report=report)
    mem_mon.insert_spikes(st_mon, value=V_th)

    return {
            'mu_offs': mu_offs,
            'mu_amp': mu_amp,
            'sigma_offs': sigma_offs,
            'sigma_amp': sigma_amp,
            'freq': freq,
            'mem': mem_mon.values,
            'spikes': st_mon.spiketimes,
            #'mu': mu_mon.values,
            #'sigma': sigma_mon.values,
            #'xi': I_mon.values,
            'duration': defaultclock.t,
            'seed': seed,
            }

if __name__=='__main__':
    data_dir = sys.argv[1]
    data = DataManager(data_dir)
    print('\n')
    mu_offs = mo_min+rand(nsims)*(mo_max-mo_min)
    mu_amp = ma_min+rand(nsims)*(ma_max-ma_min)
    sigma_offs = so_min+rand(nsims)*(so_max-so_min)
    sigma_amp = sa_min+rand(nsims)*(sa_max-sa_min)
    freq = (freq_min+rand(nsims)*(freq_max-freq_min)).round()
    params = zip(mu_offs, mu_amp, sigma_offs, sigma_amp, freq,)
    print("Simulations configured. Running ...")
    run_tasks(data, ousim, params, gui=False, poolsize=0,
            numitems=nsims)
    print("Simulations done!")

