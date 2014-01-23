from brian import *

duration = 10000*ms
dt = defaultclock.dt
spike_peak = 10*mV

gL = 1./(100*Mohm)
R = 1./gL
VL = -70*mV
Vth = -60*mV
Vr = -70*mV
C = 1e-10*farad
tau_m = C/gL

mu_amp = 0.01*namp
mu_offs = 0.1*namp
sigma_amp = 0.01*namp
sigma_offs = 0.01*namp

freq_ang = 3*Hz

eqs=Equations('dV/dt = (mu-V+VL)/tau_m + sigma*I/tau_m : volt')
eqs+=Equations('dmu/dt = (-mu_amp*sin(t*freq_ang)*freq_ang)/gL : volt')
eqs+=Equations('dI/dt = -I/dt + xi*tau_m**-0.5 : 1')
eqs+=Equations('dsigma/dt = (-sigma_amp*sin(t*freq_ang)*freq_ang)/gL : volt')
eqs.prepare()
nrn = NeuronGroup(1, eqs, reset='V=Vr', threshold='V>Vth')
nrn.mu = (mu_offs+mu_amp)/gL
nrn.sigma = (sigma_offs+sigma_amp)/gL
nrn.V = VL

mumon = StateMonitor(nrn, 'mu', record=True)
sigmamon = StateMonitor(nrn, 'sigma', record=True)

Vmon = StateMonitor(nrn, 'V', record=True)
spikemon = SpikeMonitor(nrn)
run(duration)
Vmon.insert_spikes(spikemon, spike_peak)
plot(Vmon.times, Vmon[0])
plot(Vmon.times, ones(len(Vmon[0]))*Vth)
#figure(2)
#plot(mumon.times, mumon[0], label='mu')
#plot(sigmamon.times, sigmamon[0], label='sigma')
#legend()
show()

