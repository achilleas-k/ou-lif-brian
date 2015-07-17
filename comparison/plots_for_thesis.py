import numpy as np
from brian import mV, ms, sqrt, Hz
import matplotlib.pyplot as plt

goodrmsth = 0.001
badrmsth  = 0.003

goodsdth = 0.05
badsdth  = 0.08


def get_rand(data, good):
    # rmslist = np.array([d["rms"] for d in data.itervalues()])
    sdlist  = np.array([d["sd"]  for d in data.itervalues()])
    vthlist = np.array([c[5] for c in data.iterkeys()])
    muplist = np.array([c[0]+c[1] for c in data.iterkeys()])
    muolist = np.array([c[1] for c in data.iterkeys()])
    general_rule = (vthlist == 10*mV) & (muplist >= 1.0*mV/ms) &\
        (muplist <= 2.5*mV/ms) & (muolist <= 1.5*mV/ms)
    if good:
        # valididx = np.flatnonzero((rmslist < goodrmsth) & general_rule)
        valididx = np.flatnonzero((sdlist < goodsdth) & general_rule)
    else:
        # valididx = np.flatnonzero((rmslist > badrmsth)  & general_rule)
        valididx = np.flatnonzero((sdlist > badsdth)  & general_rule)

    # fuggit --- let's plot all valids
    print("Plotting {} valid traces.".format(valididx))
    for idx in valididx:
        c = data.keys()[idx]
        print(c)
        plot_traces(data, c, "bd_")

    print("{} valid simulations.".format(len(valididx)))
    if len(valididx) == 0:
        print("No valid indices. Change your thresholds!")
    randidx = np.random.choice(valididx)
    randc = data.keys()[randidx]

    print(randc)

    ma, mo, sa, so, freq, vth = randc
    randcsp = (ma, mo, sa, so, freq, 100*mV)
    print(randcsp)

    return randc, randcsp

def make_fname(config):
    ma, mo, sa, so, freq, vth = config
    ma = float(ma)
    mo = float(mo)
    sa = float(sa*np.sqrt(0.001)*1000)
    so = float(so*np.sqrt(0.001)*1000)
    freq = int(freq)
    vth = int(vth*1000)
    basename = "ma_{}-mo_{}-sa_{}-so_{}-f_{}-vth_{}".format(
        ma, mo, sa, so, freq, vth).replace(".", "_")
    return basename+".pdf"

def plot_traces(data, config, prefix):
    d   = data[config]
    t   = d["OU"]["t"]
    ou  = d["OU"]["V"]
    lif = d["LIF"]["V"]
    plt.figure()
    plt.plot(t, ou*1000, color="b")
    plt.plot(t, lif*1000, color="g")
    plt.plot(t, np.zeros_like(t)+10, "k--")
    plt.axis([0, 0.5, 0, 20])
    # nsp = len(d["OU"]["spikes"])+len(d["LIF"]["spikes"])
    plt.xlabel("t (s)")
    plt.ylabel("V (mV)")
    plt.savefig(prefix+make_fname(config))

data = np.load("results.npz")["data"].item()

# confignsp, configsp = get_rand(data, good=True)
# plot_traces(data, confignsp, "gd_")
# plot_traces(data, configsp,  "gd_")
# confignsp, configsp = get_rand(data, good=False)
# plot_traces(data, confignsp, "bd_")
# plot_traces(data, configsp,  "bd_")

# nice good config
config = (0.5*mV/ms, 1.0*mV/ms, 0.1*mV/sqrt(ms), 0.1*mV/sqrt(ms), 10*Hz, 100*mV)
plot_traces(data, config, "gd_")
config = (0.5*mV/ms, 1.0*mV/ms, 0.1*mV/sqrt(ms), 0.1*mV/sqrt(ms), 10*Hz, 10*mV)
plot_traces(data, config, "gd_")

print("Stats for 'good' sample:")
d = data[config]
print("Spike distance: {}\nMaximum difference: {}\nRMS: {}".format(
    d["sd"], d["md"], d["rms"]))

config = (0.5*mV/ms, 1.0*mV/ms, 0.1*mV/sqrt(ms), 0.5*mV/sqrt(ms), 20*Hz, 100*mV)
plot_traces(data, config, "bd_")
config = (0.5*mV/ms, 1.0*mV/ms, 0.1*mV/sqrt(ms), 0.5*mV/sqrt(ms), 20*Hz, 10*mV)
plot_traces(data, config, "bd_")

print("Stats for 'bad' sample:")
d = data[config]
print("Spike distance: {}\nMaximum difference: {}\nRMS: {}".format(
    d["sd"], d["md"], d["rms"]))
