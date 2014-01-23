import brian_no_units
from brian import *
from brian.library.random_processes import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import itertools
import neurotools as nt
import sys
import os
import gc

duration = 2*second
N = 30 # number of simulations per configuration
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
    sigma_offs = sigma_offs*mV
    sigma_amp = sigma_amp*mV
    freq = freq*Hz
    freq_ang = freq*2*pi # angluar frequency for equation
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
    #mem_mon = StateMonitor(group, 'V', record=True)
    st_mon = SpikeMonitor(group)
    run(duration, report=None)
    #mem_mon.insert_spikes(st_mon, value=v_th)
    #mslope = 0
    #slopes = array([])
    #for i in range(N):
    #    cmslope, cslopes = nt.firing_slope(mem_mon[i], st_mon[i])
    #    mslope += cmslope
    #    slopes = append(slopes, cslopes)
    #mslope /= N
    spiketimes = array([])
    for n in xrange(N):
        spiketimes = append(spiketimes, st_mon[n])
    return {
            'mu_offs': mu_offs,
            'mu_amp': mu_amp,
            'sigma_offs': sigma_offs,
            'sigma_amp': sigma_amp,
            'freq': freq,
            #'mem': mem_mon.values,
            'spikes': spiketimes,
            #'mslope': mslope,
            #'slopes': slopes
            }

def hist_spike_density(data, yvarname):
    yvars = array([])
    spikes_per = array([])
    miny = sys.maxint
    maxy = 0
    for d in data.itervalues():
        spikes = d.get('spikes')
        period = 1./d.get('freq')
        yvar = d.get(yvarname)
        yvar *= 1000
        yvars = append(yvars, ones(len(spikes))*yvar)
        spikes_per = append(spikes_per, mod(spikes, period))
        miny = yvar if yvar < miny else miny
        maxy = yvar if yvar > maxy else maxy
    maxy = 40   # Y max value override
    heatmap, xedges, yedges = histogram2d(yvars, spikes_per,
            bins=[linspace(miny, maxy, 20), linspace(0, period, 120)])
    extent = [0, period, miny, maxy]
    '''
    Using vmax works as a cuttoff point for the maximum spike density
    to show. Everything above is considered to be max.
    It is like cutting off the peaks of a 3D plot and regarding everything
    above the cut height to be maximum or above.
    '''
    imshow(heatmap, aspect='auto', origin='lower', extent=extent,
            vmin=0, vmax=20, cmap=mpl.cm.gray_r)

def scatter_spikes(data, yvarname):
    yvars = array([])
    spikes_per = array([])
    miny = sys.maxint
    maxy = 0
    for d in data.itervalues():
        spikes = d.get('spikes')
        period = 1./d.get('freq')
        yvar = array(d.get(yvarname))
        yvar *= 1000
        yvars = ones(len(spikes))*yvar
        spikes_per = mod(spikes, period)
        scatter(spikes_per, yvars, c='b')
        miny = yvar if yvar < miny else miny
        maxy = yvar if yvar > maxy else maxy
    maxy = 40   # Y max value override
    axis([0, period, miny, maxy])

if __name__=='__main__':
    no_sigma_data_dir = "stochres_no_sigma_mu_amp"
    const_sigma_data_dir = "stochres_const_sigma_mu_amp"
    sin_sigma_data_dir = "stochres_sin_sigma_mu_amp"
    mu_offs = [0]
    mu_amp = linspace(0, 60, 61)
    sigma_offs = [5]
    sigma_amp = [5]
    freq = [5]
    variable = 'mu_amp'
    config_string = "mo%i_ma%i_so%i_sa%i_f%i" % (max(mu_offs), max(mu_amp),
            max(sigma_offs), max(sigma_amp), max(freq))
    print "Created individual iterators"
    params_no_sigma = itertools.product(mu_offs, mu_amp,
            [0], [0], freq)
    params_const_sigma = itertools.product(mu_offs, mu_amp,
            add(sigma_offs, sigma_amp), [0], freq)
    params_sin_sigma = itertools.product(mu_offs, mu_amp,
            sigma_offs, sigma_amp, freq)
    print "Created configuration generator (product)"
    no_sigma_data = DataManager(no_sigma_data_dir)
    const_sigma_data = DataManager(const_sigma_data_dir)
    sin_sigma_data = DataManager(sin_sigma_data_dir)
    print("Simulations configured. Running ...")
    run_tasks(no_sigma_data, ousim, params_no_sigma, gui=False)
    run_tasks(const_sigma_data, ousim, params_const_sigma, gui=False)
    run_tasks(sin_sigma_data, ousim, params_sin_sigma, gui=False)
    print("Simulations done. Generating plots.")

    f1 = figure(1)
    f1.clf()
    suptitle(r'$\mu_o$ = 10 mV, $\mu_A$ = 1-40 mV, f = 5 Hz')
    subplot(411)
    title("No noise")
    hist_spike_density(no_sigma_data, variable)
    subplot(412)
    title(r'Constant noise, $\sigma$ = 10 mV')
    hist_spike_density(const_sigma_data, variable)
    subplot(413)
    title(r'Sinusoidal noise, $\sigma_o$ = 5 mV, $\sigma_A$ = 5 mV')
    hist_spike_density(sin_sigma_data, variable)
    period = 1./freq[0]
    subplot(414)
    title("Signal")
    t = arange(0, period+0.01, 0.01)
    freq_ang = freq[0]*Hz*2*pi
    plot(t,sin(t*freq_ang))
    yticks([])
    axis([0,period,-1,1])
    xlabel("t (sec)")
    subplot(413)
    ylabel(r'Signal amplitude ($\mu$, mV)', {'verticalalignment' : 'bottom'})
    subplots_adjust(hspace=0.6, right=0.85)
    colorax = f1.add_axes((0.9, 0.1, 0.05, 0.8))
    colorbar(cax=colorax)
    savefig("figures/"+config_string+"_heatmap.png")

    f2 = figure(2)
    f2.clf()
    suptitle(r'$\mu_o$ = 10 mV, $\mu_A$ = 1-40 mV, f = 5 Hz')
    subplot(411)
    title("No noise")
    scatter_spikes(no_sigma_data, variable)
    subplot(412)
    title(r'Constant noise, $\sigma$ = 10 mV')
    scatter_spikes(const_sigma_data, variable)
    subplot(413)
    title(r'Sinusoidal noise, $\sigma_o$ = 5 mV, $\sigma_A$ = 5 mV')
    scatter_spikes(sin_sigma_data, variable)
    period = 1./freq[0]
    subplot(414)
    title("Signal")
    t = arange(0, period+0.01, 0.01)
    freq_ang = freq[0]*Hz*2*pi
    plot(t,sin(t*freq_ang))
    yticks([])
    axis([0,period,-1,1])
    xlabel("t (sec)")
    subplot(413)
    ylabel(r'Signal amplitude ($\mu$, mV)', {'verticalalignment' : 'bottom'})
    subplots_adjust(hspace=0.6)
    savefig("figures/"+config_string+"_scatter.png")
