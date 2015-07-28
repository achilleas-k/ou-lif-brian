'''
Plot all the figures for the 2013 Brain Research paper as they will appear
in the manuscript.
'''
from __future__ import division
import matplotlib.pyplot as plt
from matplotlib.mlab import frange
# mpl.use('Agg')
import numpy as np

#rc('font',**{'family':'sans-serif','sans-serif':['Bitstream Vera Sans']})
#rc('text', usetex=True) # use latex
plot_save_loc = '/home/achilleas/tmp'

def plot_est_vs_act_err(x_data, y_data, x_label, y_label, figname='',
                        axis_limits=None, plot_zeros=False, size=(800,600),
                        plot_title='', ticks=None):
    '''
    Used to draw the seven parameter estimation plots:
    - Frequency
    - mu_p, mu_o, mu_a
    - sigma_p, sigma_o, sigma_a
    '''
    fontsize = 20
    if not figname:
        figname = '%s_vs_%s' % (x_label, y_label)
    if not len(x_data) or not len(y_data):
        return
    if not plot_zeros:
        nz_i = np.nonzero(x_data)
        y_data = np.array(y_data)[nz_i]
        x_data = np.array(x_data)[nz_i]

    '''
    Compute averages and error bars.
    '''
    x_data = x_data.round(4) # rounding prevents duplicates due to fp precision
    y_data = y_data.round(4)
    positions = []
    data = []
    # organise data for box plots
    for ux in np.unique(x_data):
        inds = np.flatnonzero(x_data == ux)
        positions.append(ux)
        data.append(y_data[inds])

    if figname=="frequency_estimation":
        plt.hist(np.abs(x_data-y_data), log=True)
        plt.show()

    oneline = [min(min(y_data), min(x_data)), max(max(y_data), max(x_data))]
    dpi = 100
    size = (size[0]/dpi, size[1]/dpi)
    plt.figure(figsize=size, dpi=dpi)
    plt.boxplot(data, positions=positions, widths=0.1, showfliers=False)
    plt.plot(oneline, oneline, color='gray', linestyle='--', label='1:1')
    plt.xlabel(x_label, fontsize=fontsize)
    if ticks is None:
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
    else:
        plt.xticks(ticks[0], ticks[1], fontsize=fontsize)
        plt.yticks(ticks[0], ticks[1], fontsize=fontsize)
    plt.ylabel(y_label, fontsize=fontsize)
    # title(plot_title, fontsize=fontsize)
    leftright = 0.1*12/size[0]
    plt.subplots_adjust(left=leftright, bottom=0.2, right=1-leftright)
    if axis_limits is not None:
        plt.axis(axis_limits)
    savetypes = ["pdf"] #, "png", "eps", "svg"]
    for ft in savetypes:
        filename = "%s/%s.%s" % (plot_save_loc, figname, ft)
        plt.savefig(filename)
        print("Figures saved: %s" % (filename))
    plt.clf()
    return


if __name__=='__main__':
    # Load the data
    filename = '/media/olddisks/simdata/ou_lif/lotsofdata.npz'
    archive = np.load(filename)
    mu_offs_est = archive['mu_offs_est']
    mu_offs_actual = archive['mu_offs_actual']
    mu_amp_est = archive['mu_amp_est']
    mu_amp_actual = archive['mu_amp_actual']
    mu_peaks_est = archive['mu_peaks_est']
    mu_peaks_actual = archive['mu_peaks_actual']
    mu_t_est = archive['mu_t_est']
    mu_rng_est = archive['mu_rng_est']
    sigma_offs_est = archive['sigma_offs_est']
    sigma_offs_actual = archive['sigma_offs_actual']
    sigma_amp_est = archive['sigma_amp_est']
    sigma_amp_actual = archive['sigma_amp_actual']
    sigma_peaks_est = archive['sigma_peaks_est']
    sigma_peaks_actual = archive['sigma_peaks_actual']
    freq_est = archive['freq_est']
    freq_actual = archive['freq_actual']
    spikecounts = archive['spikecounts']
    CV = archive['CV']
    mu_peaks_est_bad = archive['mu_peaks_est_bad']
    outrates = archive['outrates']
    spiketimes = archive['spiketimes']

    # Start plotting
    plot_est_vs_act_err(freq_actual, freq_est, '$f \mathrm{(Hz)}$',
                        '$\hat{f} \mathrm{(Hz)}$', 'frequency_estimation',
                        plot_zeros=True, axis_limits=[0, 25, 0, 25])

    plot_est_vs_act_err(mu_peaks_actual, mu_peaks_est,
                        '$\mu_p \mathrm{(mV/ms)}$',
                        '$\hat{\mu}_p \mathrm{(mV/ms)}$',
                        'mu_peaks_estimation',
                        plot_zeros=True,
                        axis_limits=[0, 4.5, 0, 4.5], size=(1200, 450),
                        plot_title='A'+' '*110,
                        ticks=(frange(0, 4.5, 0.5),
                               [0, '', 1, '', 2, '', 3, '', 4, '']),
                        )

    plot_est_vs_act_err(mu_offs_actual, mu_offs_est,
                        '$\mu_0 \mathrm{(mV/ms)}$',
                        '$\hat{\mu}_0 \mathrm{(mV/ms)}$', 'mu_offs_estimation',
                        plot_zeros=True,
                        axis_limits=[0, 2.5, 0, 2.5],
                        size=(600, 450),
                        plot_title='B'+' '*65,
                        ticks=(frange(0, 2.5, 0.25),
                               [0, '', 0.5, '', 1, '', 1.5, '', 2]),
                        )

    plot_est_vs_act_err(mu_amp_actual, mu_amp_est, '$\mu_a \mathrm{(mV/ms)}$',
                        '$\hat{\mu}_a \mathrm{(mV/ms)}$',
                        'mu_amp_estimation', plot_zeros=True,
                        axis_limits=[0, 2.5, 0, 2.5],
                        size=(600, 450),
                        plot_title='C'+' '*65,
                        ticks=(frange(0, 2.5, 0.25),
                               [0, '', 0.5, '', 1, '', 1.5, '', 2]),
                        )

    plot_est_vs_act_err(sigma_peaks_actual, sigma_peaks_est,
                        '$\sigma_p \mathrm{(mV/\sqrt{ms})}$',
                        '$\hat{\sigma}_p \mathrm{(mV/\sqrt{ms})}$',
                        'sigma_peaks_estimation',
                        plot_zeros=True, axis_limits=[-0.2, 2.5, -0.2, 2.5],
                        size=(1200, 450), plot_title='A'+' '*110,
                        ticks=(frange(0, 2.5, 0.25),
                               [0, '', 0.5, '', 1, '', 1.5, '', 2]),
                        )

    plot_est_vs_act_err(sigma_offs_actual, sigma_offs_est,
                        '$\sigma_0 \mathrm{(mV/\sqrt{ms})}$',
                        '$\hat{\sigma}_0 \mathrm{(mV/\sqrt{ms})}$',
                        'sigma_offs_estimation',
                        plot_zeros=True, axis_limits=[-0.2, 1.2, -0.2, 1.2],
                        size=(600, 450),
                        plot_title='B'+' '*65,
                        ticks=(frange(0, 1, 0.25), [0, '', 0.5, '', 1]),
                        )

    plot_est_vs_act_err(sigma_amp_actual, sigma_amp_est,
                        '$\sigma_a \mathrm{(mV/\sqrt{ms})}$',
                        '$\hat{\sigma}_a \mathrm{(mV/\sqrt{ms})}$',
                        'sigma_amp_estimation',
                        plot_zeros=True, axis_limits=[-0.2, 1.2, -0.2, 1.2],
                        size=(600, 450),
                        plot_title='C'+' '*65,
                        ticks=(frange(0, 1, 0.25), [0, '', 0.5, '', 1]),
                        )
