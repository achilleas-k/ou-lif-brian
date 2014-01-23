'''
Contruct 5 spike trains as follows:
    - Pure poisson spike train
    - Regular bursts at regular intervals
    - Regular bursts at irregular intervals
    - Irregular bursts at regular intervals
    - Irregular bursts at irregular intervals

Check the power spectrum of each.
Determine how to best discover the underlying frequency of the inter-burst
times, which should correspond to the signal frequency.
'''

from brian import *
import random

def find_burst_length(spiketrain):
    '''
    Estimate burst length as average number of consecutive spikes with an
    inter-spike interval < mean ISI'''
    ISIs = diff(spiketrain)
    mISI = mean(ISIs)
    sub_mISI = flatnonzero(ISIs < mISI)
    cons = 0
    ncons = []
    for ii in diff(sub_mISI):
        if ii == 1:
            cons += 1
        else:
            ncons.append(cons)
            cons = 0
    return int(mean(ncons)+1)

def plot_psd_with_components(psd, inter, intra, blength, bin, plot_title):
    burst_interval = inter+intra*blength
    inter_component = len(psd)*bin/burst_interval
    intra_component = len(psd)*bin/intra
    plot(psd, color='blue', label='PSD', zorder=5)
    plot([inter_component]*2, [min(psd), max(psd)], color='red',
            linestyle='--', label='Inter-burst frequency', linewidth=2)
    plot([intra_component]*2, [min(psd), max(psd)], color='green',
            linestyle='--', label='Intra-burst frequency', linewidth=2)
    axis(xmax=intra_component*1.8)
    legend(loc='best')
    title(plot_title)
    savefig(plot_title+'.png')
    clf()
    return

inter = 0.070*second
intra = 0.005*second
burst_length = 6
regb_regi_int = []
regb_iregi_int = []
iregb_regi_int = []
iregb_iregi_int = []
for ii in range(20):
    regb_regi_int.append(inter)
    regb_iregi_int.append(inter+random.gauss(0, float(inter)/4)*second)
    iregb_regi_int.append(inter)
    iregb_iregi_int.append(inter+random.gauss(0, float(inter)/4)*second)
    for jj in range(burst_length):
        regb_regi_int.append(intra)
        regb_iregi_int.append(intra)
        iregb_regi_int.append(intra+random.gauss(0, float(intra)/4)*second)
        iregb_iregi_int.append(intra+random.gauss(0, float(intra)/4)*second)

regb_regi_int = array(regb_regi_int)
regb_regi_st = cumsum(regb_regi_int)
regb_iregi_int = array(regb_iregi_int)
regb_iregi_st = cumsum(regb_iregi_int)
iregb_regi_int = array(iregb_regi_int)
iregb_regi_st = cumsum(iregb_regi_int)
iregb_iregi_int = array(iregb_iregi_int)
iregb_iregi_st = cumsum(iregb_iregi_int)


lamb = 1./mean(regb_regi_int)
poiss_int = []
for ii in range(len(regb_regi_int)):
    poiss_int.append(random.expovariate(lamb))
poiss_int = array(poiss_int)
poiss_st = cumsum(poiss_int)

width = 10*second
bin = 1*ms
regb_regi_acorr = autocorrelogram(regb_regi_st, width=width, bin=bin)
iregb_regi_acorr = autocorrelogram(iregb_regi_st, width=width, bin=bin)
regb_iregi_acorr = autocorrelogram(regb_iregi_st, width=width, bin=bin)
iregb_iregi_acorr = autocorrelogram(iregb_iregi_st, width=width, bin=bin)
poiss_acorr = autocorrelogram(poiss_st, width=width, bin=bin)
regb_regi_acorr = regb_regi_acorr[argmax(regb_regi_acorr):]
regb_iregi_acorr = regb_iregi_acorr[argmax(regb_iregi_acorr):]
iregb_regi_acorr = iregb_regi_acorr[argmax(iregb_regi_acorr):]
iregb_iregi_acorr = iregb_iregi_acorr[argmax(iregb_iregi_acorr):]
poiss_acorr = poiss_acorr[argmax(poiss_acorr):]
regb_regi_psd = fft(regb_regi_acorr)
regb_iregi_psd = fft(regb_iregi_acorr)
iregb_regi_psd = fft(iregb_regi_acorr)
iregb_iregi_psd = fft(iregb_iregi_acorr)
poiss_psd = fft(poiss_acorr)

regb_regi_blength = find_burst_length(regb_regi_st)
regb_iregi_blength = find_burst_length(regb_iregi_st)
iregb_regi_blength = find_burst_length(iregb_regi_st)
iregb_iregi_blength = find_burst_length(iregb_iregi_st)
poiss_blength = find_burst_length(poiss_st)

from IPython import embed
embed()

#plot_psd_with_components(regb_regi_psd, inter, intra, burst_length, bin,
#        'Regular burst at Regular intervals')
#plot_psd_with_components(regb_iregi_psd, inter, intra, burst_length, bin,
#        'Regular burst at Irregular intervals')
#plot_psd_with_components(iregb_regi_psd, inter, intra, burst_length, bin,
#        'Irregular burst at Regular intervals')
#plot_psd_with_components(iregb_iregi_psd, inter, intra, burst_length, bin,
#        'Irregular burst at Irregular intervals')
#
#plot_psd_with_components(poiss_psd, inter, intra, burst_length, bin,
#        'Poisson')

