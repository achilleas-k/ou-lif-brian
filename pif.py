from brian import *

duration = 200*second
N_sims = 1

lif_eq = ['dV/dt = I/tau_m : volt']
V_rest = 0*mV
V_reset = 0*mV
V_th = 15*mV
t_refr = 2*ms
tau_m = 10*msecond
N_in = 10
f_in = 10*Hz
DV_s = 0.3*mV
I = 0.08*mV
nrns = NeuronGroup(N_sims,lif_eq,threshold=V_th,reset=V_reset,\
        refractory=t_refr)

mem = StateMonitor(nrns, 'V', record=True)
st = SpikeMonitor(nrns)
run(duration, report='stdout')

plot(mem.times, mem[0], mem.times, ones(len(mem.times))*V_th)
title('Membrane')
show()


