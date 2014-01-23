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
w = 2*ms

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
    v_th = random()*30*mV
    tau = random()*30*ms
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
            'v_th' : v_th,
            'tau' : tau,
            'spikes': st_mon.spiketimes,
            'mslope': mslope,
            'slopes': slopes
            }


if __name__=='__main__':
    runsims = True
    clargs = sys.argv
    if len(clargs) > 1:
        if '--no-run' in clargs:
            clargs.remove('--no-run')
            runsims = False
        else:
            runsims = True
        if '--del-old' in clargs:
            clargs.remove('--del-old')
            for fname in data.session_filenames():
                os.remove(fname)
            if not runsims:
                print 'Deleted old data and set to no run.\nNo data available.'
                sys.exit()
        if len(clargs) < 2:
            print 'Must supply at least a data directory name.'
            sys.exit(1)
        else:
            data_dir = sys.argv[-1]
    else:
        print 'Must supply at least a data directory name.'
        sys.exit(1)
    if runsims:
        data_basepath = data_dir+'.data'
        if os.path.isdir(data_basepath):
            while not (quit == 'q' or quit == 'a' or quit == ''):
                quit = raw_input("Data directory already exists [%s].\
Append data or quit [a/Q]: " % os.path.abspath(data_basepath))
                quit = quit.lower()
            if quit.lower() == 'q' or quit == '':
                sys.exit(0)
        data = DataManager(data_dir)
        num_confs = 2000
        mo_min = 10; mo_max = 50
        mu_offs = random(num_confs)*(mo_max-mo_min)+mo_min
        ma_min = 10; ma_max = 50
        mu_amp = random(num_confs)*(ma_max-ma_min)+ma_min
        so_min = 0; so_max = 50
        sigma_offs = random(num_confs)*(so_max-so_min)+so_min
        sa_min = 0; sa_max = 50
        sigma_amp = random(num_confs)*(sa_max-sa_min)+sa_min
        freq_min = 10; freq_max = 10
        freq = random(num_confs)*(freq_max-freq_min)+freq_min
        configs = zip(mu_offs, mu_amp, sigma_offs, sigma_amp, freq)
        print("Simulations configured. Data will be saved to %s\nRunning ..." %\
                os.path.abspath(data.basepath))
        run_tasks(data, ousim, configs, gui=False, verbose=False)
        print("Simulations done!")
    else:
        data = DataManager(data_dir)
        print 'Data will be loaded from %s' % os.path.abspath(data.basepath)

    print "mo, ma, so, sa, fr, mslope, mslope_w, stdslope, stdslope_w, misi"
    for result in data.itervalues():
        mo = result.get('mu_offs')
        ma = result.get('mu_amp')
        so = result.get('sigma_offs')
        sa = result.get('sigma_amp')
        fr = result.get('freq')
        #tau = result.get('tau')
        #v_th = result.get('v_th')
        mslope = result.get('mslope')
        slopes = result.get('slopes')
        spikes = result.get('spikes').get(0) # assumes 1 sim per config
        if len(spikes) > 2:
            misi = mean(diff(spikes))
        else:
            misi = 0

        '''
        Calculate windowed slope value.
        '''
        mem = result.get('mem')[0]
        spikes = result.get('spikes')[0]
        mslope_w, slopes_w = nt.firing_slope(
                mem,
                spikes,
                dt=dt,
                w=w)

        '''
        Change to print only when ma and sa are zero to check data for
        classic OU LIF
        '''
        print '%f, '*10 % (mo, ma, so, sa, fr, mslope, mslope_w, std(slopes), std(slopes_w), misi)

