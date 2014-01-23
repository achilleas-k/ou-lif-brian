import brian_no_units
from brian import *
from brian.library.random_processes import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import itertools
import neurotools as nt
import sys
import os
from datetime import datetime

'''
Global variables for plot limits and alignment.
'''
xmax = 90
xmin = -5
hspace = 0.4
wspace = 0.4
top = 0.9
bottom = 0.1
left = 0.1
right = 0.9
dt = 0.1*ms

def plot_slope_isi_vs(param_name, param_arr, mslope, slopes, spikes):
    title('Slope and ISI Vs %s' % param_name)
    xlabel(param_name + " (mV)")
    mslopes_param = zeros(0)
    stdslopes_param = zeros(0)
    misis_param = zeros(0)
    stdisis_param = zeros(0)
    mcv_param = zeros(0)
    stdcv_param = zeros(0)

    sort_ind = argsort(param_arr)
    mslope = mslope[sort_ind]
    spikes = spikes[sort_ind]
    param_arr = param_arr[sort_ind]
    slopes = slopes[sort_ind]
    xvalues = unique(param_arr)
    for p in xvalues:
        '''
        There are two ways to consider calculating the standard deviation for
        the slope and the ISI.
        One is to calculate the standard deviation using all individual values.
        The other is to calculate it using only the means of each simulation.
        In the latter case, we would have the stanrd deviation of the mean of
        each quantity. It sounds weird, but it would make sense, since we
        assume that the mean slope and ISI are more indicative of what the
        neuron is doing overall, in the long-duration limit.
        '''
        p_ind = where(param_arr == p)
        mslopes_p = mean(mslope[p_ind])
        '''
        stdslopes_p = std(mslope[p_ind])
        '''
        slopes_p = []
        for s in slopes[p_ind]:
            slopes_p.extend(s)
        stdslopes_p = std(slopes_p)
        spiketrains_p = spikes[p_ind]   # array of dictionaries
        isis_p = []
        cvs_p = []
        for st_dict in spiketrains_p:
            for st in st_dict.values():
                if (len(st) > 1):
                    isis_p.append(mean(diff(st)))
                    cvs_p.append(std(diff(st))/mean(diff(st)))
                else:
                    isis_p.append(0)
                    cvs_p.append(0)
                isis_p.extend(diff(st))
        misis_p = mean(isis_p)
        stdisis_p = std(isis_p)
        mcv_p = mean(cvs_p)
        stdcv_p = std(cvs_p)
        '''Append all calculated quantities to respective arrays'''
        mslopes_param = append(mslopes_param, mslopes_p)
        stdslopes_param = append(stdslopes_param, stdslopes_p)
        misis_param = append(misis_param, misis_p)
        stdisis_param = append(stdisis_param, stdisis_p)
        mcv_param = append(mcv_param, mcv_p)
        stdcv_param = append(stdcv_param, stdcv_p)

    '''Slope line'''
    errorbar(xvalues*1000, mslopes_param, yerr=stdslopes_param, marker='.',
            linestyle='-', markersize=15, label='Mean slope')
    ylabel('Firing slope (mV/ms)')
    legend(loc=2)

    #'''ISI line'''
    #twinx()
    #errorbar(xvalues*1000, misis_param, yerr=stdisis_param, marker='^',
    #        linestyle='--', color='r', markersize=10, label='Mean ISI')
    #ylabel('Inter-Spike Interval (sec)')

    #'''CV line'''
    #twinx()
    #stdcv_param = 0
    #errorbar(xvalues*1000, mcv_param, yerr=stdcv_param, marker='^',
    #    linestyle='--', color='r', markersize=10, label='CV(ISI)')
    #ylabel('CV(ISI)')

    '''Second line legend'''
    legend(loc=1)

    '''Plot ranges'''
    axis(xmin=xmin, xmax=xmax)

if __name__=='__main__':
    if len(sys.argv) > 1:
        data_dir = sys.argv[1]
    else:
        print "Gimme a data directory name (without the .data extension)."
        sys.exit(1)
    data = DataManager(data_dir)
    '''
    Make individual plots for each input parameter vs slope and vs firing.
    Take mean slope, mean firing freq., stdev slope, stdev firing frequency.
    Newer simulations save data as list of dictionaries.
    Dictionary keys correspond to variable names.
    '''
    mu_offs, mu_amp, sigma_offs, sigma_amp, freq,\
    mslope, slopes, spikes, mem =\
            array(zip(*[[d.get('mu_offs'), d.get('mu_amp'),
                d.get('sigma_offs'), d.get('sigma_amp'),
                d.get('freq'), d.get('mslope')/dt,
                d.get('slopes')/dt, d.get('spikes'), d.get('mem')]
                for d in data.itervalues()]))

    print "Data loaded. Performing extra calculations ... "
    mslope_w = mslope.copy()
    slopes_w = slopes.copy()
    for i in range(len(mem)):
        mslope_w[i], slopes_w[i] = nt.norm_firing_slope(
                mem[i][0],
                spikes[i][0],
                th=20*mV,
                tau=20*ms,
                dt=dt,
                w=2.0*ms)

    print "Done processin. Plotting figures ..."
    figure()
    subplot(3, 2, 1)
    plot_slope_isi_vs('Mu offs', mu_offs, mslope_w, slopes_w, spikes)
    subplot(3, 2, 2)
    plot_slope_isi_vs('Mu amp', mu_amp, mslope_w, slopes_w, spikes)
    subplot(3, 2, 3)
    plot_slope_isi_vs('Sigma offs', sigma_offs, mslope_w, slopes_w, spikes)
    subplot(3, 2, 4)
    plot_slope_isi_vs('Sigma amp', sigma_amp, mslope_w, slopes_w, spikes)
    subplot(3, 2, 5)
    sigma_peak = add(sigma_offs, sigma_amp)
    plot_slope_isi_vs('Sigma peak', sigma_peak, mslope_w, slopes_w, spikes)
    subplot(3, 2, 6)
    mu_peak = add(mu_offs, mu_amp)
    plot_slope_isi_vs('Mu peak', mu_peak, mslope_w, slopes_w, spikes)
    #figure("Frequency")
    #plot_slope_isi_vs('Frequency', freq, mslope, slopes, spikes)
    subplots_adjust(left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace)
    show()
