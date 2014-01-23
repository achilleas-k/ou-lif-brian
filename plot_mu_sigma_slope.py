from brian import *
from brian.tools.datamanager import *
import itertools
import sys


if len(sys.argv) == 1:
    print "Gimme data name"
    sys.exit(1)

dataname = sys.argv[1]
data = DataManager(dataname)

numitems = data.itemcount()
mu = zeros(numitems)
sigma = zeros(numitems)
slope = zeros(numitems)
fout = zeros(numitems)

i = 0
for d in data.itervalues():
    mu_amp = d.get('mu_amp')
    mu_offs = d.get('mu_offs')
    sigma_amp = d.get('sigma_amp')
    sigma_offs = d.get('sigma_offs')
    #spikes = d.get('spikes')
    #ISIs = array([])
    #for sp in spikes.itervalues():
    #    ISIs = append(ISIs, diff(sp))

    mu_peak = mu_amp+mu_offs
    sigma_peak = sigma_amp+sigma_offs

    mslope = d.get('mslope')

    mu[i] = mu_peak
    sigma[i] = sigma_peak
    slope[i] = mslope
    #fout[i] = 1./mean(ISIs) if len(ISIs) > 1 else 0
    i += 1
    sys.stdout.write("\r%4.1f %% complete" % (100*float(i)/numitems))
    sys.stdout.flush()

mu_unique = unique(mu)
sigma_unique = unique(sigma)

slope_means = zeros([len(mu_unique), len(sigma_unique)])
fout_means = zeros([len(mu_unique), len(sigma_unique)])
slope_stds = zeros([len(mu_unique), len(sigma_unique)])
fout_stds = zeros([len(mu_unique), len(sigma_unique)])

for y, x in itertools.product(mu_unique, sigma_unique):
    ypos = where(mu_unique == y)[0]
    xpos = where(sigma_unique == x)[0]

    inds = logical_and(mu == y, sigma == x)

    slope_means[ypos, xpos] = mean(slope[inds])
    #fout_means[ypos, xpos] = mean(fout[inds])

    slope_stds[ypos, xpos] = std(slope[inds])
    #fout_stds[ypos, xpos] = std(fout[inds])




extent = [sigma_unique[0], sigma_unique[-1], mu_unique[0], mu_unique[-1]]

figure(1)
title("Mean slope")
imshow(slope_means, extent=extent, origin='lower', interpolation='gaussian')
xlabel("sigma")
ylabel("mu")
colorbar()
#figure(2)
#title("Stdev slope")
#imshow(slope_stds, extent=extent, origin='lower', interpolation='gaussian')
#colorbar()
#xlabel("sigma")
#ylabel("mu")

#figure(3)
#title("Mean firing rate")
#imshow(fout_means, extent=extent, origin='lower', interpolation='gaussian')
#colorbar()
#xlabel("sigma")
#ylabel("mu")
#figure(4)
#title("Stdev firing rate")
#imshow(fout_stds, extent=extent, origin='lower', interpolation='gaussian')
#colorbar()
#xlabel("sigma")
#ylabel("mu")
show()
