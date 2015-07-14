import numpy as np
import matplotlib.pyplot as plt

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

if __name__ == "__main__":
    data = np.load("results.npz")["data"].item()
    maxsd = find_max_item(data, "sd")
    maxmd = find_max_item(data, "md")
    maxsq = find_max_item(data, "sq")
    plot_item(maxsd, "maxsd")
    plot_item(maxmd, "maxmd")
    plot_item(maxsq, "maxsq")
