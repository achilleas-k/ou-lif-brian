from brian import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import sys
import gc
from time import time
import neurotools as nt
import itertools

forever = 1e10*second
N = 10
defaultclock.dt = dt = 0.1*ms
V_th = 10*mV
tau = 10*ms
t_refr = 0*ms
v_reset = 0*mV
V0 = 0*mV
min_spikes = 5000

def ousim(mu_offs, mu_amp, sigma_offs, sigma_amp, freq, report):
    clear(True)
    gc.collect()
    reinit_default_clock()
    np.random.seed(int(time()+(mu_offs+mu_amp+sigma_offs+sigma_amp)*1e3))
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
    st_mon = SpikeMonitor(group)
    mem_mon = StateMonitor(group, 'V', record=True)

    duration = 10*second
    #duration = forever
    #@network_operation
    #def stop_condition():
    #    if st_mon.nspikes >= min_spikes*1.0/N:
    #        stop()

    run(duration, report=report)
    mem_mon.insert_spikes(st_mon, V_th)
    return {
            'mu_offs': mu_offs,
            'mu_amp': mu_amp,
            'sigma_offs': sigma_offs,
            'sigma_amp': sigma_amp,
            'freq': freq,
            'spikes': st_mon.spiketimes,
            #'mem': mem_mon,
            'duration': defaultclock.t,
            }

def psd_freq_est(st):
    '''
    Estimate the frequency based on the fourier transform of the
    autocorrelogram (power spectrum) of the spike train
    '''
    width = 1
    bin = 1*ms
    acorr = autocorrelogram(st, width=width, bin=bin, T=width)
    acorr = acorr[len(acorr)/2:] # acorr is symmetric
    psd = fft(acorr)
    cutfreq = len(psd)/2
    psd[cutfreq:] = 0 # low pass filter
    psd[0] = 0 # remove first component
    max_amp_i = argmax(psd)
    max_freq = (1.0*max_amp_i/len(psd))/bin
    return max_freq, psd

if __name__=='__main__':
    data_dir = 'noisesearch'
    data = DataManager(data_dir)
    print('\n')
    if '--no-run' not in sys.argv:
        mu_offs = [0.6]
        mu_amp = [0.5]
        sigma_offs = frange(0.0, 0.1, 0.01)
        sigma_amp = [0]
        freq = [5, 10, 20]
        params_prod = itertools.product(mu_offs,
                                        mu_amp,
                                        sigma_offs,
                                        sigma_amp,
                                        freq,)
        nsims = len(mu_offs)*len(mu_amp)*\
                len(sigma_offs)*len(sigma_amp)*\
                len(freq)

        print("Simulations configured. Running ...")
        run_tasks(data, ousim, params_prod, gui=False, poolsize=0,
                numitems=nsims)
        print("Simulations done!")
    else:
        print("Skipping simulation run. Working with %s.data directory\n" % (
            data_dir))

    numsims = data.itemcount()
    print("Total number of simulations: %i" % numsims)
    print("Estimating sine wave frequencies ... ")

    mu_offs = []
    mu_amp = []
    sigma_offs = []
    sigma_amp = []
    freq = []
    freq_est = []
    freq_est_err = []
    all_psd = []
    spikecounts = []
    for ind, d in enumerate(data.itervalues()):
        i = ind+1
        sys.stdout.write("\rProcessing %i of %i" % (i, numsims))
        sys.stdout.flush()
        spikes = d.get('spikes')
        spikes = d.get('spikes')
        mo = d.get('mu_offs')
        ma = d.get('mu_amp')
        so = d.get('sigma_offs')
        sa = d.get('sigma_amp')
        fr = d.get('freq')
        duration = d.get('duration')

        for st in spikes.values():
            freq_est_cur = 1e42*Hz
            psd = []
            try:
                freq_est_cur, psd = psd_freq_est(st)
            except Exception:
                print "\n\nException was raised at %i with f: %s, so: %s" % (
                        ind, fr, so)
            freq_est.append(freq_est_cur)
            error = abs(freq_est_cur-fr)/fr
            mu_offs.append(mo)
            mu_amp.append(ma)
            sigma_offs.append(so)
            sigma_amp.append(sa)
            freq.append(fr)
            freq_est_err.append(error)
            spikecounts.append(len(st))
            #all_psd.append(psd)

    archivename = 'ns_results.npz'
    np.savez(archivename,
        mu_offs=mu_offs,
        mu_amp=mu_amp,
        sigma_offs=sigma_offs,
        sigma_amp=sigma_amp,
        freq=freq,
        freq_est=freq_est,
        freq_est_err=freq_est_err,
        spikecounts=spikecounts,
        )
    print("\nResults have been saved in the %s archive" % archivename)
    sys.exit(1)


#################################################################
    print("Drawing PSD plots")
    sorted_ind = sorted(range(len(sigma_offs)), key=sigma_offs.__getitem__)
    figure(figsize=(16, 12), dpi=100)
    for ii, si in enumerate(sorted_ind):
        sys.stdout.write('.')
        sys.stdout.flush()
        subplot(4, 4, (ii+1) % 16)
        psd = all_psd[si]
        fr_act = freq[si]
        fr_act_ind = int(fr_act*len(psd)*ms)
        fr_est = freq_est[si]
        fr_est_ind = int(fr_est*len(psd)*ms)
        plot(psd)
        plot([fr_act_ind]*2, [min(psd), max(psd)], color='green',
                linestyle='--', label='Actual fr')
        plot([fr_est_ind]*2, [min(psd), max(psd)], color='red',
                linestyle='--', label='Estim. fr')
        so =  sigma_offs[si]*sqrt(ms)/mV
        title('so: %2.2f, act: %2.0f, est: %2.0f' % (so, fr_act, fr_est))
        if (ii+1) % 16 == 0:
            savefig('psds_%i.png' % ii)
            sys.stdout.write(' figure saved\n')
            clf()


'''
Use the following to load the data in another script or an interactive
session:

filename = 'ns_results.npz'
archive = np.load(filename)
mu_offs = archive['mu_offs']
mu_amp = archive['mu_amp']
sigma_offs = archive['sigma_offs']
sigma_amp = archive['sigma_amp']
freq = archive['freq']
freq_est = archive['freq_est']
freq_est_err = archive['freq_est_err']
spikecounts = archive['spikecounts']

'''
