from brian import *

tau = 100*ms
xthreshold = inf*mV
dt = defaultclock.dt
muv = -15*mV
mux = -15*mV
seegma = 10*mV
eqs=Equations('dy/dt=-y/dt + xi*tau**-0.5 : 1')
eqs+=Equations('dx/dt=(mux-x)/tau+sigma*y/dt : volt')
eqs+=Equations('dsigma/dt = cos(t*10*Hz)*200*mV/dt : volt')
eqs+=Equations('dv/dt = (-v+muv)/tau + seegma*xi*tau**-0.5 : volt')
group = NeuronGroup(1000, eqs)
group.sigma = 10*mV
group.x=100*mV
group.v=100*mV

smon = StateMonitor(group, 'sigma', record=True)
xmon = StateMonitor(group, 'x', record=True)
vmon = StateMonitor(group, 'v', record=True)
run(1*second, report='stdout')

xv = xmon.values
vv = vmon.values

print "mean(x): %f, var(x): %f" % (mean(xv), var(xv))
print "mean(v): %f, var(v): %f" % (mean(vv), var(vv))
#T = size(xmon[0])
#segs = 10
#stds = [std(xmon[0][s*T/segs:(s+1)*T/segs]) for s in range(segs)]
#print stds

subplot(211)
plot(xmon.times, xmon[0])
plot(smon.times, smon[0])
title("X")
subplot(212)
plot(vmon.times, vmon[0])
title("V")
show()
