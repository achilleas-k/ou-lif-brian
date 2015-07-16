import numpy as np
import matplotlib.pyplot as plt
from brian import display_in_unit, mV, ms, Hz, sqrt

duration = 0.5  # seconds

param_names = ["$\mu_a$", "$\mu_0$", "$\sigma_a$", "$\sigma_0$",
               "$f$", "$V_{th}$"]
param_units = [mV/ms, mV/ms, mV/sqrt(ms), mV/sqrt(ms),
               Hz, mV]

def find_max_item(data, key):
    keyvals = []
    for d in data.itervalues():
        keyvals.append(d[key])
    print("{}: {}".format(key, max(keyvals)))
    return data.values()[np.argmax(keyvals)]

def plot_item(item, name):
    plt.figure(name)
    plt.plot(item["LIF"]["t"]*1000, item["LIF"]["V"]*1000)
    plt.plot(item["OU"]["t"]*1000,  item["OU"]["V"]*1000)
    plt.xlabel("t (ms)")
    plt.ylabel("V (mV)")
    plt.savefig(name+".pdf")

def colour_hist(data, cidx):
    pn = param_names[cidx]
    un = param_units[cidx]
    grouped_sd = {}
    grouped_md = {}
    grouped_sq = {}
    maxsd = 0
    maxmd = 0
    maxsq = 0
    for config in data:
        k = config[cidx]
        if k in grouped_sd:
            grouped_sd[k].append(data[config]["sd"])
            grouped_md[k].append(data[config]["md"])
            grouped_sq[k].append(data[config]["sq"])
        else:
            grouped_sd[k] = [data[config]["sd"]]
            grouped_md[k] = [data[config]["md"]]
            grouped_sq[k] = [data[config]["sq"]]
        maxsd = max(maxsd, data[config]["sd"])
        maxmd = max(maxsd, data[config]["md"])
        maxsq = max(maxsd, data[config]["sq"])

    for k in sorted(grouped_sd.iterkeys()):
        plt.figure("Spike distance histogram")
        y, x = np.histogram(grouped_sd[k], bins=np.linspace(0, maxsd, 21))
        plt.plot(x[:-1], y, label="{} = {}".format(pn, display_in_unit(k, un)))
        plt.figure("Max difference histogram")
        y, x = np.histogram(grouped_md[k], bins=np.linspace(0, maxmd, 21))
        plt.plot(x[:-1], y, label="{} = {}".format(pn, display_in_unit(k, un)))
        plt.figure("Squared diff   histogram")
        y, x = np.histogram(grouped_sq[k], bins=np.linspace(0, maxsq, 21))
        plt.plot(x[:-1], y, label="{} = {}".format(pn, display_in_unit(k, un)))

    plt.figure("Spike distance histogram")
    plt.legend()
    plt.savefig("sdhist_{}.pdf".format(pn.replace("$","")))
    plt.clf()

    plt.figure("Max difference histogram")
    plt.legend()
    plt.savefig("mdhist_{}.pdf".format(pn.replace("$","")))
    plt.clf()

    plt.figure("Squared diff   histogram")
    plt.legend()
    plt.savefig("sqhist_{}.pdf".format(pn.replace("$","")))
    plt.clf()

def plot_all_hist(data):
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
    plt.hist(square_difference[nsp]*1000/duration, bins=50)
    plt.axis(xmin=0)
    plt.xlabel("Mean square difference (mV$^2$/s)")
    plt.savefig("square_difference.pdf")


if __name__ == "__main__":
    data = np.load("results.npz")["data"].item()
    # maxsd = find_max_item(data, "sd")
    # maxmd = find_max_item(data, "md")
    # maxsq = find_max_item(data, "sq")
    # plot_item(maxsd, "maxsd")
    # plot_item(maxmd, "maxmd")
    # plot_item(maxsq, "maxsq")
    # plot_all_hist(data)

    configlen = len(data.keys()[0])
    print("Ploting colour histograms")
    for c in range(configlen):
        colour_hist(data, c)
        print("Finished with {}".format(param_names[c]))
