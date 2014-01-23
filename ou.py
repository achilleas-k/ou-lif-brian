#!/usr/bin/env python

from brian import *
from brian.library.random_processes import *
import neurotools
import sys

tau_m = 10*ms
in_sigma = float(sys.argv[2])*mA
in_mu = float(sys.argv[1])*mA
v_th = 15*mV
C = 1*mfarad
t_refr = 2*ms
v_reset = 0*mV

eqs=Equations('dv/dt=-v/tau_m+I/C : volt')
eqs+=OrnsteinUhlenbeck('I',mu=in_mu,sigma=in_sigma,tau=tau_m)
group = NeuronGroup(1, eqs, threshold=v_th, refractory=t_refr, reset=v_reset)
group.rest()

mem = StateMonitor(group, 'v', record=True)
I_mon = StateMonitor(group, 'I', record=True)
st = SpikeMonitor(group)
run(1000*ms)

(sta_avg, sta_std, sta_wins) = neurotools.sta(mem[0],st[0],2*ms,0.1*ms)
print "Spikes:",st.nspikes
print "---"
print "mu:",in_mu
print "sigma:",in_sigma
print "---"
print "Slope:",mean(diff(sta_avg))
#subplot(211)
plot(mem.times,mem[0])#,mem.times,ones(len(mem.times))*v_th)
#title('Membrane')
#subplot(212)
#plot(sta_avg)
#title('Reverse correlation (membrane)')
show()


