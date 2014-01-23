from __future__ import division
from brian import *
import sys
import gc
from time import time
import neurotools as nt

duration = 0.3*second
defaultclock.dt = dt = 0.1*ms
h = dt
V_th = 10*mV
tau = 10*ms
b = 1/tau
t_refr = 0*ms
v_reset = 0*mV
V0 = 0*mV

mu_estimates = []
mu_actuals = []

sigma_estimates = []
sigma_actuals = []

ion()
while True:
    nrepeats = raw_input('How many? ')
    nrepeats = int(nrepeats)
    for i in range(nrepeats):
        sys.stdout.write('\r%i/%i ... ' % (i+1, nrepeats))
        sys.stdout.flush()
        clear(True)
        gc.collect()
        reinit_default_clock()
        mu = round(rand()*3, 3)
        sigma = round(rand()*3, 3)
        mu = mu*mV/ms
        sigma = sigma*mV/sqrt(ms)

        eqs =Equations('dV/dt = mu-(V+V0)/tau + sigma*I/sqrt(dt) : volt')
        eqs+=Equations('dI/dt = -I/dt + xi/sqrt(dt) : 1')
        eqs.prepare()
        group = NeuronGroup(1, eqs, threshold='V>V_th', refractory=t_refr,
                reset=v_reset)
        group.V = V0
        mem_mon = StateMonitor(group, 'V', record=True)
        st_mon = SpikeMonitor(group)
        run(duration, report=None)
        if st_mon.nspikes < 5:
            continue
        mem_mon.insert_spikes(st_mon, value=V_th)
        spikes = st_mon[0]
        mem = mem_mon[0]

        mu_est = []
        sigma_est = []
        '''
        BEGIN ESTIMATION
        '''
        for tp, t in zip(spikes[:-1], spikes[1:]):
            t *= second
            tp *= second
            isi = t-tp
            K = round(isi/h)
            Vi = mem[tp/h+1:t/h+1]
            'mu estimation'
            me = V_th / (K/b*(1-exp(-b*h))) + b/K * sum(Vi[:-1])*volt

            'sigma estimation'
            se2 = 2/K * sum(((Vi[1:] - mu/b + (mu/b - Vi[:-1]) * exp(-b*h))**2)/(1/b*(1-exp(-2*b*h))))
            se = sqrt(se2)

            mu_est.append(me)
            sigma_est.append(se)
        '''
        END ESTIMATION
        '''
        mu_estimates.append(mean(mu_est))
        mu_actuals.append(mu)
        sigma_estimates.append(mean(sigma_est))
        sigma_actuals.append(sigma)
    sys.stdout.write('\n')
    clf()
    scatter(mu_actuals, mu_estimates, c='b')
    scatter(sigma_actuals, sigma_estimates, c='r')
    xlabel('act')
    ylabel('est')
    plot([0, 2], [0, 2], 'k--')
    show()
    print("mu corr. coef.: %f" % corrcoef(mu_actuals, mu_estimates)[0,1])
    print("sigma corr. coef.: %f" % corrcoef(sigma_actuals, sigma_estimates)[0,1])
    ask_for_inp = True
    while ask_for_inp:
        reply = raw_input("Would you like to know more? [Y/n/i] ")
        if reply in ['n', 'N']:
            ask_for_inp = False
            sys.exit(0)
        elif reply in ['i', 'I']:
            import IPython
            IPython.embed()
            ask_for_inp = True
        elif reply in ['y', 'Y']:
            ask_for_inp = False
        else:
            ask_for_inp = False
