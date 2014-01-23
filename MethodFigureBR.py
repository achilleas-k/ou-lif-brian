from brian import *

defaultclock.dt = dt = 0.1*ms

freq = 10*Hz
period = 1/freq
duration = 4*period

mu_offs = 3*mV/ms
mu_amp = 1*mV/ms
sigma_offs = 1*mV/sqrt(ms)
sigma_amp = 1*mV/sqrt(ms)

T = frange(0*second, duration, 0.1*ms)
signal = mu_amp*sin(freq*T*2*pi)+mu_offs
ncoef = sigma_amp*sin(freq*T*2*pi)+sigma_offs
noise = normal(0, 1, len(T))*ncoef

# let's get random samples from the original sine wave
# the larger the value, the higher the probability of getting a sample
probabilities = (signal-min(signal))/max(signal)
probabilities *= 0.3   # so that peaks are not p=1
T_est = []
# the following isn't the simplest way to do it, but allows us to check
# for refractory period
for ti, p in zip(T, probabilities):
    if rand() < p:
        if not T_est or ti-T_est[-1] >= 0.002:
            T_est.append(ti)
T_est = array(T_est)
est = mu_amp*sin(freq*T_est*2*pi)+mu_offs
rel_errors = random(len(est))*0.2-0.1
est += est*rel_errors

duraticks = frange(period, duration, period)
duraticklabels = ['$%i f^{-1}$' % (i+1) if i else '$ f^{-1} $'
                                        for i in range(len(duraticks))]

periodticks = [0, float(period)/2, float(period)]
periodticklabels = ['0', '$0.5 f^{-1}$', '$f^{-1}$']

muticks = [float(mu_offs-mu_amp), float(mu_offs), float(mu_offs+mu_amp)]
muticklabels = ['$\mu_0-\mu_a$', '$\mu_0$', '$\mu_p$']

periodspikes = T_est % period
binedges = linspace(0*second, period, 11)
avgest = []
for l, r in zip(binedges[:-1], binedges[1:]):
    binspike_inds = bitwise_and(periodspikes > l, periodspikes < r)
    binmu = est[binspike_inds]
    avgest.append(mean(binmu) if len(binmu) else 0)
avgest = array(avgest)

axeslimits_full = [0, float(duration),
        (mu_offs-mu_amp)*0.95, (mu_offs+mu_amp)*1.05]
axeslimits_period = [0, float(period),
        (mu_offs-mu_amp)*0.95, (mu_offs+mu_amp)*1.05]


titalign = ' '*110 # spaces to align title to the left

# start plotting
figure(1)
clf()
subplot(5, 1, 1)
plot(T, signal+noise, 'k-')
plot(T, signal, color='grey', linestyle='--')
title('A%s' % titalign)
#title('Original signal (with noise)')
xticks(duraticks, duraticklabels)
yticks(muticks, muticklabels)
axis(axeslimits_full)

subplot(5, 1, 2)
stem(T_est, ones(len(T_est))*(mu_offs+mu_amp/2), linefmt='k-', basefmt='k-',
        markerfmt='w-', bottom=mu_offs)
plot(T, signal, color='grey', linestyle='--')
title('B%s' % titalign)
#title('Fired spikes')
xticks(duraticks, duraticklabels)
yticks(muticks, muticklabels)
#yticks([])
axis(axeslimits_full)

subplot(5, 1, 3)
plot(T_est, est, 'k.')
plot(T, signal, color='grey', linestyle='--')
title('C%s' % titalign)
#title('Estimated $\mu$ values $\hat{\mu}_i$')
xticks(duraticks, duraticklabels)
yticks(muticks, muticklabels)
axis(axeslimits_full)

subplot(5,1,4)
plot(periodspikes, est, 'k.')
plot(T, signal, color='grey', linestyle='--')
title('D%s' % titalign)
#title('Estimated $\mu$ values, aligned to one period $\hat{\mu}_i$')
xticks(periodticks, periodticklabels)
yticks(muticks, muticklabels)
axis(axeslimits_period)

subplot(5,1,5)
plot(binedges[:-1]+0.5*period/10, avgest, color='black', linestyle='-',
        marker='.', markersize=10)
plot(T, signal, color='grey', linestyle='--')
title('E%s' % titalign)
#title('Binned average of estimated $\mu$ values $\hat{\mu}_b$')
xticks(periodticks, periodticklabels)
yticks(muticks, muticklabels)
axis(axeslimits_period)

subplots_adjust(hspace=0.9)


savefig('brain_research_figures/method.png')
savefig('brain_research_figures/method.pdf')
savefig('brain_research_figures/method.ps')
