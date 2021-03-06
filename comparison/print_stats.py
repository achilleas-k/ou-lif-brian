import numpy as np
import matplotlib.pyplot as plt
from brian import display_in_unit, mV, ms, Hz, sqrt

duration = 0.5  # seconds

param_names = ["$\mu_a$", "$\mu_0$", "$\sigma_a$", "$\sigma_0$",
               "$f$", "$V_{th}$"]
param_units = [mV/ms, mV/ms, mV/sqrt(ms), mV/sqrt(ms),
               Hz, mV]

def find_max_item(data, key):
    if key == "sd":
        spikes = True
    else:
        spikes = False
    maxvalue = 0
    for config in data:
        d = data[config]
        Vth = config[5]
        # nspikes = len(d["OU"]["spikes"])+len(d["LIF"]["spikes"])
        # if ((nspikes > 0) and spikes) or ((nspikes == 0) and not spikes):
        if ((Vth > 50*mV) and not spikes) or ((Vth < 50*mV) and spikes):
            if (config[2] > 0.3*mV/sqrt(ms)) or (config[3] > 0.3*mV/sqrt(ms)):
                continue
            if d[key] > maxvalue:
                maxvalue = d[key]
                maxkey = config
    print("{}: {} \t({})".format(key, maxvalue, maxkey))
    return data[maxkey]

def plot_item(item, name):
    plt.figure(name)
    plt.plot(item["LIF"]["t"]*1000, item["LIF"]["V"]*1000)
    plt.plot(item["OU"]["t"]*1000,  item["OU"]["V"]*1000)
    plt.xlabel("t (ms)")
    plt.ylabel("V (mV)")
    plt.savefig(name+".pdf")

def colour_hist(data, cidx, spikes):
    pn = param_names[cidx]
    un = param_units[cidx]
    grouped_sd = {}
    grouped_md = {}
    grouped_rms = {}
    maxsd = 0
    maxmd = 0
    maxsq = 0
    for config in data:
        k = config[cidx]
        if ((len(data[config]["OU"]["spikes"]) > 0) ^ spikes):
            continue
        if k in grouped_sd:
            grouped_sd[k].append(data[config]["sd"])
            grouped_md[k].append(data[config]["md"]*1e3)
            grouped_rms[k].append(data[config]["rms"]*1e3)
        else:
            grouped_sd[k] = [data[config]["sd"]]
            grouped_md[k] = [data[config]["md"]*1e3]
            grouped_rms[k] = [data[config]["sq"]*1e3]
        maxsd = max(maxsd, data[config]["sd"])
        maxmd = max(maxsd, data[config]["md"]*1e3)
        maxsq = max(maxsd, data[config]["sq"]*1e3)

    for k in sorted(grouped_sd.iterkeys()):
        plt.figure("Spike distance histogram")
        y, x = np.histogram(grouped_sd[k], bins=np.linspace(0, maxsd, 6), normed=True)
        plt.plot(x[:-1], y, label="{} = {}".format(pn, display_in_unit(k, un)))
        plt.figure("Max difference histogram")
        y, x = np.histogram(grouped_md[k], bins=np.linspace(0, maxmd, 6), normed=True)
        plt.plot(x[:-1], y, label="{} = {}".format(pn, display_in_unit(k, un)))
        plt.figure("RMSE histogram")
        y, x = np.histogram(grouped_rms[k], bins=np.linspace(0, maxsq, 6), normed=True)
        plt.plot(x[:-1], y, label="{} = {}".format(pn, display_in_unit(k, un)))

    if spikes:
        sx = "spks"
    else:
        sx = "nspk"
    plt.figure("Spike distance histogram")
    plt.legend()
    plt.xlabel("SPIKE-distance")
    plt.ylabel("Number of samples")
    plt.savefig("sdhist_{}_{}.pdf".format(pn.replace("$",""), sx))
    plt.clf()

    plt.figure("Max difference histogram")
    plt.legend()
    plt.xlabel("Maximum difference (mV)")
    plt.ylabel("Number of samples")
    plt.savefig("mdhist_{}_{}.pdf".format(pn.replace("$",""), sx))
    plt.clf()

    plt.figure("RMSE histogram")
    plt.legend()
    plt.xlabel("Root mean squared error (mV)")
    plt.ylabel("Number of samples")
    plt.savefig("rmshist_{}_{}.pdf".format(pn.replace("$",""), sx))
    plt.clf()

def plot_all_hist(data):
    spike_distance = np.array([d["sd"]  for d in data.itervalues()])
    max_difference = np.array([d["md"]  for d in data.itervalues()])
    rms            = np.array([d["rms"] for d in data.itervalues()])
    nspikes = np.array([len(d["OU"]["spikes"])+len(d["LIF"]["spikes"])
                        for d in data.itervalues()])
    sp = nspikes > 0
    nsp = nspikes == 0

    # mu_a = np.array([c[0] for c in data.iterkeys()])
    # mu_0 = np.array([c[1] for c in data.iterkeys()])

    sigma_a = np.array([c[2] for c in data.iterkeys()])
    sigma_0 = np.array([c[3] for c in data.iterkeys()])
    # small sigma
    ss = (sigma_a < 0.5*mV/sqrt(ms)) & (sigma_0 < 0.5*mV/sqrt(ms))

    idx = sp & ss
    plt.figure("Spike distance")
    plt.hist(spike_distance[idx], bins=50)
    plt.axis(xmin=0)
    plt.xlabel("SPIKE-distance")
    plt.savefig("spike_distance.pdf")

    idx = nsp & ss
    print("MAX V {}".format(max([max(d["OU"]["V"]) for d in np.array(data.values())[idx]])*1000))
    plt.figure("Max deviation")
    plt.hist(max_difference[idx]*1000, bins=50)
    plt.axis(xmin=0)
    plt.xlabel("Maximum deviation (mV)")
    plt.savefig("max_difference.pdf")

    plt.figure("Squared difference")
    plt.hist(rms[idx]*1000, bins=50)
    plt.axis(xmin=0)
    plt.xlabel("Root mean squared error (mV)")
    plt.savefig("rms.pdf")


if __name__ == "__main__":
    data = np.load("results.npz")["data"].item()
    maxsd = find_max_item(data, "sd")
    maxmd = find_max_item(data, "md")
    maxsq = find_max_item(data, "sq")
    maxrms = find_max_item(data, "rms")

    plot_item(maxsd, "maxsd")
    plot_item(maxmd, "maxmd")
    plot_item(maxsq, "maxsq")
    plot_item(maxrms, "maxrms")
    plot_all_hist(data)

    configlen = len(data.keys()[0])
    # print("Ploting colour histograms")
    # for c in range(configlen):
    #     colour_hist(data, c, False)
    #     colour_hist(data, c, True)
    #     print("Finished with {}".format(param_names[c]))
