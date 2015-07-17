import numpy as np
from brian import mV
import matplotlib.pyplot as plt

goodrmsth = 0.001
badrmsth  = 0.003


def get_rand(data, good):
    rmslist = np.array([d["rms"] for d in data.itervalues()])
    vthlist = np.array([c[5] for c in data.iterkeys()])
    if good:
        valididx = np.flatnonzero((rmslist < goodrmsth) & (vthlist == 100*mV))
    else:
        valididx = np.flatnonzero((rmslist > badrmsth)  & (vthlist == 100*mV))
    if len(valididx) == 0:
        print("No valid indices. Change your thresholds!")
    randidx = np.random.choice(valididx)
    randc = data.keys()[randidx]

    print(randc)

    ma, mo, sa, so, freq, vth = randc
    randcsp = (ma, mo, sa, so, freq, 15*mV)
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
    nsp = len(d["OU"]["spikes"])+len(d["LIF"]["spikes"])
    if (nsp > 0):
        plt.xlabel("t (s)")
    plt.ylabel("V (mV)")
    plt.savefig(prefix+make_fname(config))

data = np.load("results.npz")["data"].item()

confignsp, configsp = get_rand(data, good=True)
plot_traces(data, confignsp, "gd_")
plot_traces(data, configsp,  "gd_")
confignsp, configsp = get_rand(data, good=False)
plot_traces(data, confignsp, "bd_")
plot_traces(data, configsp,  "bd_")
