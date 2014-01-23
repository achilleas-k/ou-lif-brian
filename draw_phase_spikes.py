from brian import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import sys
import gc
from time import time
from warnings import warn

filename = sys.argv[1]
archive = np.load(filename)
mu_offs_actual = archive['mu_offs_actual']
mu_amp_actual = archive['mu_amp_actual']
mu_means_actual = archive['mu_means_actual']
mu_peaks_actual = archive['mu_peaks_actual']
phase_spikes = archive['phase_spikes']

if not len(unique(mu_offs_actual)) == 1:
    warn('Non unique mu offs')
mo = mu_offs_actual[0]
figure(figsize=(14.40, 9), dpi=100)
title('mu offs = %4.2f mV/ms' % mo)
imshow(phase_spikes[argsort(mu_amp_actual)], aspect='auto', origin='lower',
        interpolation='none')
colorbar()
yticks(range(len(mu_amp_actual)), sorted(mu_amp_actual))
ylabel('mu amp (mV/ms)')
xlabel('t (ms)')
savefig('phspikes_dma_mo%4.2f.png' % mo)
print('Image saved as phspikes_dma_mo%4.2f.png' % mo)
