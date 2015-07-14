from __future__ import print_function
import os
import brian
from brian import (Network, defaultclock, clock, Equations, NeuronGroup,
                   PulsePacket, SpikeGeneratorGroup, Connection,
                   PoissonGroup, sqrt,
                   StateMonitor, SpikeMonitor, mV, ms, second, Hz)
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from spikerlib.metrics import kreuz
from multiprocessing.pool import Pool
import itertools as it

sin = brian.sin
pi = brian.pi

duration = 1000*ms
tau = 10*ms
t_refr = 0*ms
V_reset = 0*mV
V0 = 0*mV

results_t = [0.3*second, 0.8*second]

dt = defaultclock.dt

def ousim(mu_amp, mu_offs, sigma_amp, sigma_offs, freq, V_th):
    # mu_amp, mu_offs, sigma_amp, sigma_offs, freq, V_th = config
    if sigma_amp > sigma_offs:
        sigma_amp = sigma_offs
    # print("Setting up OU LIF simulation...")
    ounet = Network()
    clock.reinit_default_clock()
    eqs =Equations('dV/dt = mu-(V+V0)/tau + sigma*I/sqrt(dt) : volt')
    eqs+=Equations('dI/dt = -I/dt + xi/sqrt(dt) : 1')
    eqs+=Equations('mu = mu_amp*sin(t*freq*2*pi) + mu_offs : volt/second')
    eqs+=Equations('sigma = sigma_amp*sin(t*freq*2*pi) + sigma_offs :'
                                                        ' volt/sqrt(second)')
    eqs.prepare()
    ounrn = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr,
                                                                reset=V_reset)
    ounet.add(ounrn)

    ounrn.V = V0
    V_mon = StateMonitor(ounrn, 'V', record=True)
    st_mon = SpikeMonitor(ounrn)
    ounet.add(V_mon, st_mon)

    ounet.run(duration)

    V_mon.insert_spikes(st_mon, value=V_th*2)
    times = V_mon.times
    membrane = V_mon[0]
    return times, st_mon.spiketimes[0], membrane

def lifsim(mu_amp, mu_offs, simga_amp, sigma_offs, freq, V_th):
    lifnet = Network()
    clock.reinit_default_clock()
    eqs = Equations('dV/dt = (-V+V0)/tau : volt')
    eqs.prepare()
    lifnrn = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr,
                         reset=V_reset)
    lifnet.add(lifnrn)
    pulse_times = (np.arange(1, duration*freq, 1)+0.25)/freq
    pulse_spikes = []
    Npoiss = 5000
    Npulse = 5000
    wpoiss = (mu_offs-mu_amp)/(Npoiss*freq)
    wpulse = mu_amp/(Npulse*freq)
    sigma = 1/(freq*5)
    if (wpulse != 0):
        for pt in pulse_times:
            pp = PulsePacket(t=pt*second, n=Npulse, sigma=sigma)
            pulse_spikes.extend(pp.spiketimes)
        pulse_input = SpikeGeneratorGroup(Npulse, pulse_spikes)
        pulse_conn = Connection(pulse_input, lifnrn, 'V', weight=wpulse)
        lifnet.add(pulse_input, pulse_conn)
    if (wpoiss != 0):
        poiss_input = PoissonGroup(Npoiss, freq)
        poiss_conn = Connection(poiss_input, lifnrn, 'V', weight=wpoiss)
        lifnet.add(poiss_input, poiss_conn)

    V_mon = StateMonitor(lifnrn, 'V', record=True)
    st_mon = SpikeMonitor(lifnrn)
    lifnet.add(V_mon, st_mon)

    lifnet.run(duration)

    V_mon.insert_spikes(st_mon, value=V_th*2)
    times = V_mon.times
    membrane = V_mon[0]
    return times, st_mon.spiketimes[0], membrane

def process_results(ou, lif, config):
    times_ou, spikes_ou, voltage_ou = ou
    times_lif, spikes_lif, voltage_lif = lif

    start, end = results_t
    end -= start
    start -= start
    dist = kreuz.distance(spikes_lif, spikes_ou,
                          start, end, (end-start)/(1*ms))
    kdist = np.trapz(dist[1], dist[0])
    maxdiff = np.max(np.abs(voltage_lif-voltage_ou))
    sqdiff = np.sum(np.square(voltage_lif-voltage_ou))

    data = load_data("results.npz")
    OUdict  = {"V": voltage_ou,   "t": times_ou,   "spikes": spikes_ou}
    LIFdict = {"V": voltage_lif,  "t": times_lif,  "spikes": spikes_lif}
    data[config] = {"OU": OUdict, "LIF": LIFdict,
                    "sd": kdist, "md": maxdiff, "sq": sqdiff}
    np.savez("results.npz", data=data)

def make_plots(ou, lif, fnamesuffix):
    times_ou, spikes_ou, voltage_ou = ou
    times_lif, spikes_lif, voltage_lif = lif
    start, end = results_t
    end -= start
    start -= start
    ax_limits= {"xmin": start*1000, "xmax": end*1000,
                "ymin": min(0, min(voltage_ou*1000)),
                "ymax": max(0, max(voltage_ou*1000))}

    plt.figure(figsize=(8, 6))
    plt.subplot2grid((5, 1), (0, 0), rowspan=4, colspan=1)
    plt.plot(times_lif*1000, voltage_lif*1000)
    plt.plot(times_ou*1000, voltage_ou*1000)
    plt.axis(**ax_limits)
    plt.ylabel("mV")
    xt, _ = plt.xticks()
    plt.xticks(xt, [])
    plt.subplot2grid((5, 1), (4, 0), rowspan=1, colspan=1)
    plt.plot(times_ou*1000, abs(voltage_ou-voltage_lif)*1000)
    plt.xlabel("t (ms)")
    plt.ylabel("mV")
    plt.yticks([0, 2.5, 5])
    plt.axis(**ax_limits)
    plt.axis(ymax=5)
    mpl.rcParams["font.size"] = 12
    plt.subplots_adjust(left=0.1, top=0.95, bottom=0.1, right=0.95, hspace=0.2)
    plt.savefig("ou_vs_lif_"+fnamesuffix+".pdf")

    plt.figure(figsize=(8, 3))
    plt.plot(times_ou*1000, voltage_ou*1000)
    plt.xlabel("t (ms)")
    plt.ylabel("mV")
    plt.axis(**ax_limits)
    mpl.rcParams["font.size"] = 12
    plt.subplots_adjust(left=0.1, top=0.95, bottom=0.2, right=0.95)
    plt.savefig("ou_sin_"+fnamesuffix+".pdf")

    plt.figure(figsize=(8, 3))
    plt.plot(times_lif*1000, voltage_lif*1000)
    plt.xlabel("t (ms)")
    plt.ylabel("mV")
    plt.axis(**ax_limits)

    mpl.rcParams["font.size"] = 12
    plt.subplots_adjust(left=0.1, top=0.95, bottom=0.2, right=0.95)
    plt.savefig("lif_sin_"+fnamesuffix+".pdf")

def cut_results(results, start, end):
    time, spikes, voltage = results
    tidx = np.flatnonzero((start <= time) & (time < end))
    time = time[tidx]
    voltage = voltage[tidx]
    spikes = spikes[(start <= spikes) & (spikes < end)]
    time -= start
    spikes -= start
    return time, spikes, voltage

def load_data(fname):
    if os.path.exists(fname):
        return np.load(fname)["data"].item()
    else:
        return {}

if __name__=='__main__':
    pool = Pool()
    # param values
    mu_amp     = [0.5*mV/ms, 1.0*mV/ms, 1.5*mV/ms, 2.0*mV/ms]
    mu_offs    = [0.5*mV/ms, 1.0*mV/ms, 1.5*mV/ms, 2.0*mV/ms]
    sigma_amp  = [0.1*mV/sqrt(ms)]#, 0.5*mV/sqrt(ms), 1.0*mV/sqrt(ms)]
    sigma_offs = [0.1*mV/sqrt(ms), 0.5*mV/sqrt(ms)]#, 1.0*mV/sqrt(ms)]
    freq       = [10*Hz, 20*Hz]
    V_th       = [10*mV, 15*mV, 100*mV]

    configs = it.product(mu_amp, mu_offs, sigma_amp, sigma_offs, freq, V_th)
    data = load_data("results.npz")
    configs = [c for c in configs if c not in data]
    print("{} configurations total".format(len(configs)))
    # using apply_async instead of map to see progress
    pool_ou =  [pool.apply_async(ousim,  c) for c in configs]
    pool_lif = [pool.apply_async(lifsim, c) for c in configs]
    results_ou = []
    results_lif = []
    print("Running OU...")
    for res in pool_ou:
        res.wait()
        cres = res.get()
        cres = cut_results(cres, *results_t)
        results_ou.append(cres)
        print("{} of {} complete".format(len(results_ou), len(pool_ou)))
    print("Running LIF...")
    for res in pool_lif:
        res.wait()
        cres = res.get()
        cres = cut_results(cres, *results_t)
        results_lif.append(cres)
        print("{} of {} complete".format(len(results_lif), len(pool_lif)))

    print("Processing simulation data...")
    for idx in range(len(configs)):
        process_results(results_ou[idx], results_lif[idx], configs[idx])
        print("{} of {} complete".format(idx+1, len(configs)))
    data = load_data("results.npz")
    spike_distance    = np.array([d["sd"] for d in data.itervalues()])
    max_difference    = np.array([d["md"] for d in data.itervalues()])
    square_difference = np.array([d["sq"] for d in data.itervalues()])
    nspikes = np.array([len(d["OU"]["spikes"])+len(d["LIF"]["spikes"])
                        for d in data.itervalues()])
    sp = nspikes > 0
    nsp = nspikes == 0

    plt.figure("Spike distance")
    plt.hist(spike_distance[sp], bins=50)
    plt.axis(xmin=0)
    plt.xlabel("SPIKE-distance")
    plt.savefig("spike_distance.pdf")

    plt.figure("Max deviation")
    plt.hist(max_difference[nsp]*1000, bins=50)
    plt.axis(xmin=0)
    plt.xlabel("Maximum deviation (mV)")
    plt.savefig("max_difference.pdf")

    plt.figure("Squared difference")
    plt.hist(square_difference[nsp]*1000, bins=50)
    plt.axis(xmin=0)
    plt.xlabel("Summed square difference (mV$^2$)")
    plt.savefig("square_difference.pdf")
