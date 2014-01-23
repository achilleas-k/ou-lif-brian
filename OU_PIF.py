import brian_no_units_no_warnings
from brian import *
from brian.library.random_processes import *
import neurotools
import sys

duration = 10*second
tau_m = 10*ms
in_sigma = float(sys.argv[2])*mA
in_mu = float(sys.argv[1])*mA
V_th = 15*mV
C = 1*mfarad
t_refr = 2*ms
V_reset = 0*mV

eqs=Equations('dV/dt = in_mu/tau_m + in_sigma*xi : 1')
group = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr, reset=V_reset)

mem_mon = StateMonitor(group, 'V', record=True)
st_mon = SpikeMonitor(group)
run(duration, report='stdout')

(sta_avg, sta_std, sta_wins) = neurotools.sta(mem_mon[0],st_mon[0],0.2*ms)
print "Spikes:",st_mon.nspikes
print "---"
print "mu:",in_mu
print "sigma:",in_sigma
print "---"
print "Slope:",mean(diff(sta_avg))
plot(mem_mon.times,mem_mon[0],mem_mon.times,ones(len(mem_mon.times))*V_th)
title('Membrane')
show()


