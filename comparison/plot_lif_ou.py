from __future__ import print_function
from brian import (Network, defaultclock, clock, Equations, NeuronGroup,
                   PulsePacket, SpikeGeneratorGroup, Connection,
                   PoissonGroup, display_in_unit,
                   StateMonitor, SpikeMonitor, mV, ms, second, Hz)
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi
import matplotlib as mpl
from spikerlib.metrics import kreuz
from multiprocessing.pool import Pool

duration = 1000*ms
tau = 20*ms
t_refr = 0*ms
V_reset = 0*mV
V0 = 0*mV

mu_amp = 2.2*mV/ms
mu_offs = 2.2*mV/ms
sigma_amp = 0.1*mV/sqrt(ms)
sigma_offs = 0.1*mV/sqrt(ms)
freq = 10*Hz
freq_ang = freq*2*pi

results_t = [0.3*second, 0.8*second]

dt = defaultclock.dt

def ousim(V_th):
    print("Setting up OU LIF simulation...")
    ounet = Network()
    clock.reinit_default_clock()
    eqs =Equations('dV/dt = mu-(V+V0)/tau + sigma*I/sqrt(dt) : volt')
    eqs+=Equations('dI/dt = -I/dt + xi/sqrt(dt) : 1')
    eqs+=Equations('mu = mu_amp*sin(t*freq_ang) + mu_offs : volt/second')
    eqs+=Equations('sigma = sigma_amp*sin(t*freq_ang) + sigma_offs :'
                                                        ' volt/sqrt(second)')
    eqs.prepare()
    ounrn = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr,
                                                                reset=V_reset)
    ounet.add(ounrn)

    ounrn.V = V0
    print("mu_amp:    {} mV/ms,       mu_offs:    {} mV/ms\n"
          "sigma_amp: {} mV/sqrt(ms), sigma_offs: {} mV/sqrt(ms)\n"
          "frequency: {} Hz".format(
              mu_amp, mu_offs, sigma_amp*sqrt(ms)/mV, sigma_offs*sqrt(ms)/mV,
              freq))
    print("Configuring monitors...")
    V_mon = StateMonitor(ounrn, 'V', record=True)
    st_mon = SpikeMonitor(ounrn)
    ounet.add(V_mon, st_mon)

    print("Running simulation for {} seconds".format(duration))
    ounet.run(duration)
    print("OU done")

    V_mon.insert_spikes(st_mon, value=V_th*2)
    times = V_mon.times
    membrane = V_mon[0]
    return times, st_mon.spiketimes[0], membrane

def lifsim(V_th):
    print("Setting up LIF simulation...")
    lifnet = Network()
    clock.reinit_default_clock()
    eqs = Equations('dV/dt = (-V+V0)/tau : volt')
    eqs.prepare()
    lifnrn = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr,
                                                                reset=V_reset)
    lifnet.add(lifnrn)
    pulse_times = (np.arange(1, duration*freq, 1)+0.25)/freq
    pulse_spikes = []
    print("Generating input spike trains...")
    Nin = 5000
    mu_peak = mu_offs+mu_amp
    weight = mu_offs/Nin/freq
    Npoiss = int(Nin*(mu_offs-mu_amp)/mu_peak)
    Npulse = Nin-Npoiss
    sigma = 1/(freq*5)
    if (Npulse > 0):
        for pt in pulse_times:
            pp = PulsePacket(t=pt*second, n=Npulse, sigma=sigma)
            pulse_spikes.extend(pp.spiketimes)
        pulse_input = SpikeGeneratorGroup(Npulse, pulse_spikes)
        pulse_conn = Connection(pulse_input, lifnrn, 'V', weight=weight)
        lifnet.add(pulse_input, pulse_conn)
    if (Npoiss > 0):
        poiss_input = PoissonGroup(Npoiss, freq)
        poiss_conn = Connection(poiss_input, lifnrn, 'V', weight=weight)
        lifnet.add(poiss_input, poiss_conn)

    print("N_in: {}, Npulses: {}\n"
          "S_in: {}, sigma_in: {} ms, weight: {} mV".format(
              Nin, len(pulse_times), 1.0, sigma*1000, weight*1000))

    print("Configuring monitors...")
    V_mon = StateMonitor(lifnrn, 'V', record=True)
    st_mon = SpikeMonitor(lifnrn)
    lifnet.add(V_mon, st_mon)

    print("Running simulation for {} seconds".format(duration))
    lifnet.run(duration)
    print("LIF done")

    V_mon.insert_spikes(st_mon, value=V_th*2)
    times = V_mon.times
    membrane = V_mon[0]
    return times, st_mon.spiketimes[0], membrane

def process_results(ou, lif):
    times_ou, spikes_ou, voltage_ou = ou
    times_lif, spikes_lif, voltage_lif = lif

    start, end = results_t
    end -= start
    start -= start
    dist = kreuz.distance(spikes_lif, spikes_ou,
                          start, end, (end-start)/(1*ms))
    kdist = np.trapz(dist[1], dist[0])
    print("Spike train distance  : {}".format(kdist))
    maxdiff = np.max(np.abs(voltage_lif-voltage_ou))
    print("Max mem potential diff: {}".format(maxdiff))
    sqdiff = np.sum(np.square(voltage_lif-voltage_ou))
    print("Sum sq potential diff : {}".format(sqdiff))

def make_plots(ou, lif, fnamesuffix):
    times_ou, spikes_ou, voltage_ou = ou
    times_lif, spikes_lif, voltage_lif = lif
    start, end = results_t
    end -= start
    start -= start
    ax_limits= [start*1000, end*1000, 0*mV, (mu_amp+mu_offs)*tau*1100]

    plt.figure(figsize=(8, 6))
    plt.subplot2grid((5, 1), (0, 0), rowspan=4, colspan=1)
    plt.plot(times_lif*1000, voltage_lif*1000)
    plt.plot(times_ou*1000, voltage_ou*1000)
    plt.axis(ax_limits)
    plt.ylabel("mV")
    xt, _ = plt.xticks()
    plt.xticks(xt, [])
    plt.subplot2grid((5, 1), (4, 0), rowspan=1, colspan=1)
    plt.plot(times_ou*1000, abs(voltage_ou-voltage_lif)*1000)
    plt.xlabel("t (ms)")
    plt.ylabel("mV")
    plt.yticks([0, 2.5, 5])
    plt.axis(ax_limits)
    plt.axis(ymax=5)
    mpl.rcParams["font.size"] = 12
    plt.subplots_adjust(left=0.1, top=0.95, bottom=0.1, right=0.95, hspace=0.2)
    plt.savefig("ou_vs_lif_"+fnamesuffix+".pdf")

    plt.figure(figsize=(8, 3))
    plt.plot(times_ou*1000, voltage_ou*1000)
    plt.xlabel("t (ms)")
    plt.ylabel("mV")
    plt.axis(ax_limits)
    mpl.rcParams["font.size"] = 12
    plt.subplots_adjust(left=0.1, top=0.95, bottom=0.2, right=0.95)
    plt.savefig("ou_sin_"+fnamesuffix+".pdf")

    plt.figure(figsize=(8, 3))
    plt.plot(times_lif*1000, voltage_lif*1000)
    plt.xlabel("t (ms)")
    plt.ylabel("mV")
    plt.axis(ax_limits)

    mpl.rcParams["font.size"] = 12
    plt.subplots_adjust(left=0.1, top=0.95, bottom=0.2, right=0.95)
    plt.savefig("lif_sin_"+fnamesuffix+".pdf")

def cut_results(results, start, end):
    time, spikes, voltage = results
    tidx = np.flatnonzero((start <= time) & (time < end))
    time = time[tidx]
    voltage = voltage[tidx]
    spikes = spikes[(start <= spikes) & (spikes < end)]
    # make start == 0
    time -= start
    spikes -= start
    return time, spikes, voltage


if __name__=='__main__':
    pool = Pool()
    V_th = [5*mV, 10*mV, 15*mV, 100*mV]
    V_th = [100*mV]
    poolres = []
    poolres.append(pool.map_async(ousim, V_th))
    poolres.append(pool.map_async(lifsim, V_th))
    results = []
    for res in poolres:
        res.wait()
        cres = res.get()
        cres = [cut_results(cr, *results_t) for cr in cres]
        results.extend(cres)

    for idx in range(len(V_th)):
        suffix = display_in_unit(V_th[idx], mV).replace(" ", "_")
        process_results(results[idx], results[idx+len(V_th)])
