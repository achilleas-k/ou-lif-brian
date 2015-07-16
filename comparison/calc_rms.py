import numpy as np

data = np.load("results.npz")["data"].item()

for config in data:
    simulation = data[config]
    ouv  =  simulation["OU"]["V"]
    lifv =  simulation["LIF"]["V"]

    rms = np.sqrt(np.mean((ouv-lifv)**2))
    simulation["rms"] = rms
    data[config] = simulation
np.savez("results.npz", data=data)
