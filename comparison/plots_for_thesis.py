import numpy as np
from brian import mV
import matplotlib.pyplot as plt

goodrmsth = 0.0015
badrmsth  = 0.0035


def get_rand(data, good):
    # good here means RMS < 1.5 mV
    configs = data.keys()
    ridx = np.random.choice(len(configs))
    randc = configs[ridx]
    if good:
        while not (data[randc]["rms"] < goodrmsth and randc[5] == 100*mV):
            ridx = np.random.choice(len(configs))
            randc = configs[ridx]
    else:
        while not (data[randc]["rms"] > badrmsth  and randc[5] == 100*mV):
            ridx = np.random.choice(len(configs))
            randc = configs[ridx]

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

def plot_traces(data, config):
    d   = data[config]
    t   = d["OU"]["t"]
    ou  = d["OU"]["V"]
    lif = d["LIF"]["V"]
    if (data[config]["rms"] > badrmsth):
        prefix = "bd_"
    else:
        prefix = "gd_"
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
plot_traces(data, confignsp)
plot_traces(data, configsp)
confignsp, configsp = get_rand(data, good=False)
plot_traces(data, confignsp)
plot_traces(data, configsp)
