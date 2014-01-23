from brian import *
from brian.tools.datamanager import *
from neurotools import times_to_bin
import sys
import os
import gc
'''
Find best averaging window by taking moving averages of the spike train using
a range of values and using them to estimate the underlying frequency.
Minimising the estimation error and observing the best window values should
provide some insight.


---
There's definitely some sort of memory leak here.
'''
data = DataManager('montecarlo_10us')
dt = float(0.1*ms)
t = frange(0, 3, dt)
freq_act = []
freq_est = []
best_wins = []
for d in data.itervalues():
    spikes = d.get('spikes')
    mu_offs = d.get('mu_offs')
    mu_amp = d.get('mu_amp')
    sigma_offs = d.get('sigma_offs')
    sigma_amp = d.get('sigma_amp')
    freq = d.get('freq')
    filehead = "tmp/m%f-%f_s%f-%f-f%f" % (mu_offs/(mV/ms),
                    mu_amp/(mV/ms),
                    sigma_offs/(mV/sqrt(ms)),
                    sigma_amp/(mV/sqrt(ms)),
                    freq/Hz)
    sys.stdout.write('='*100+'\n')
    sys.stderr.write('='*100+'\n')

    for st in spikes.itervalues():
        if len(st) < 5:
            continue
        binst = times_to_bin(st)
        isis = diff(st)
        winscales = linspace(1*ms, 100*ms, 100)/dt
        all_win_estimates = []
        for ws in winscales:
            sys.stderr.write('.')
            sys.stderr.flush()
            x = frange(-ws*4, ws*4, 1) # 4 stds per side
            win = normpdf(x, 0, ws)
            frate = convolve(binst, win, mode='valid')
            fft_rate = fft(frate)
            fft0 = append(0, fft_rate[1:])
            max_amp_i = argmax(fft0)
            freq_est_cur = max_amp_i/(len(fft0)*dt)
            all_win_estimates.append(freq_est_cur)
        all_win_errors = abs(subtract(all_win_estimates,freq))
        best_est_i = argmin(all_win_errors)
        best_est = all_win_estimates[best_est_i]
        best_win = winscales[best_est_i]
        freq_est.append(best_est)
        freq_act.append(freq)
        best_wins.append(best_win)
        printstring = '\rAct: %f, Est: %f [%f, %f], Win: %f [%f, %f]\n' % (
            freq, best_est, min(all_win_estimates), max(all_win_estimates),
            best_win, min(winscales), max(winscales))
        sys.stdout.write(printstring)
        sys.stdout.flush()
        sys.stderr.write(printstring)
        clear(True)
        gc.collect()
