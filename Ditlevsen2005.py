#import brian_no_units
from brian import *
from brian.library.random_processes import *
from brian.tools.taskfarm import *
from brian.tools.datamanager import *
import itertools
import neurotools as nt
import sys
import os
from time import time
from datetime import datetime
import gc

duration = 1e10*second
N = 10 # number of simulations per configuration-process
dt = defaultclock.dt
tau = 10*ms
v_th = 10*mV
t_refr = 0*ms
v_reset = 0*mV
V0 = 0*mV
min_spikes = 1000

def ousim(mu_offs, mu_amp, sigma_offs, sigma_amp, freq, report):
    clear(True)
    gc.collect()
    reinit_default_clock()
    np.random.seed(int(time()+(mu_offs+mu_amp+sigma_offs+sigma_amp)*1e10))
    mu_offs = mu_offs*mvolt/ms
    mu_amp = mu_amp*mvolt/ms
    sigma_amp = sigma_amp*mvolt/sqrt(ms)
    sigma_offs = sigma_offs*mvolt/sqrt(ms)\
            if sigma_offs*mvolt/sqrt(ms) > sigma_amp else sigma_amp
    freq = freq*Hz
    freq_ang = freq*2*pi # angluar frequency for equation

    mu = mu_offs
    sigma = sigma_offs
    eqs=Equations('dV/dt = mu-(V+V0)/tau + sigma*xi : volt')
    #eqs+=Equations('dI/dt = -I/dt + xi*dt**-0.5 : 1')
    #eqs+=Equations('mu = mu_amp*sin(t*freq_ang)+mu_offs : volt/second')
    #eqs+=Equations('sigma = sigma_amp*sin(t*freq_ang)+sigma_offs : volt/second')
    eqs.prepare()
    group = NeuronGroup(N, eqs, threshold=v_th, refractory=t_refr,
            reset=v_reset)
    group.V = V0
    #mu_mon = StateMonitor(group, 'mu', record=True)
    #sigma_mon = StateMonitor(group, 'sigma', record=True)
    mem_mon = StateMonitor(group, 'V', record=True)
    st_mon = SpikeMonitor(group)

    @network_operation
    def stop_condition():
        if st_mon.nspikes/N >= min_spikes:
            stop()

    run(duration, report=report)
    #subplot(311)
    #plot(mu_mon.times, mu_mon[0])
    #subplot(312)
    #plot(sigma_mon.times, sigma_mon[0])
    #mem_mon.insert_spikes(st_mon, value=v_th)
    #mslope = 0
    #slopes = array([])
    #for i in range(N):
    #    cmslope, cslopes = nt.firing_slope(mem_mon[i], st_mon[i], w=2*ms)
    #    mslope += cmslope
    #    slopes = append(slopes, cslopes)
    #mslope /= N
    print "Sim with mu %s and sigma %s fired %i spikes in %s" % (
            display_in_unit(mu_offs, mV/ms),
            display_in_unit(sigma_offs, mV/sqrt(ms)),
            st_mon.nspikes,
            defaultclock.t)
    return {
            'mu_offs': mu_offs,
            'mu_amp': mu_amp,
            'sigma_offs': sigma_offs,
            'sigma_amp': sigma_amp,
            'freq': freq,
            'mem': mem_mon,
            'spikes': st_mon,
            #'mslope': mslope,
            #'slopes': slopes
            }


if __name__=='__main__':
    data_dir = sys.argv[1]
    data = DataManager(data_dir)
    print "\n\n"
    if '--no-run' not in sys.argv:
        mu_amp = [0]
        mu_offs = linspace(1, 4, 10)
        sigma_amp = [0]
        sigma_offs = linspace(0.1, 3, 10)
        freq = [10]
        nsims=len(mu_amp)*len(mu_offs)*len(sigma_amp)*len(sigma_offs)*len(freq)
        print "Created individual iterators"
        params_prod = itertools.product(mu_offs, mu_amp, sigma_offs, sigma_amp,
                freq)
        print "Created configuration generator (product)"
        print("Simulations configured. Running ...")
        run_tasks(data, ousim, params_prod, gui=True, poolsize=0,
                numitems=nsims)
        print("Simulations done!")
    else:
        print "Skipping simulation run. Working with %s.data directory" %\
                data_dir

    numsims = data.itemcount()
    i = 1
    print "Total number of simulations: %i" % numsims
    mu_est = array([])
    mu_actual = array([])
    sigma_est = array([])
    sigma_actual = array([])
    theta_est = array([])
    theta_actual = array([])

    for d in data.itervalues():
        print "\nProcessing %i of %i" % (i, numsims)
        mu_offs = d.get('mu_offs')
        mu_amp = d.get('mu_amp')
        sigma_offs = d.get('sigma_offs')
        sigma_amp = d.get('sigma_amp')
        mu = mu_offs
        sigma = sigma_offs
        filehead = "tmp/m%fs%f" % (mu, sigma)

        ss_depol = mu*tau
        if (ss_depol < v_th+0.1*mV and ss_depol > v_th-0.1*mV):
            regime = 0
            regime_str = "threshold"
        elif (ss_depol > v_th):
            regime = 1
            regime_str = "super threshold"
        elif (ss_depol < v_th):
            regime = -1
            regime_str = "sub threshold"

        print "Steady state depolarisation: %f ( %s )" % (
                ss_depol, regime_str)


        spike_mon = d.get('spikes')
        mem_mon = d.get('mem')
        #duration = mem_mon.times[-1]*second
        progress = 0
        sys.stdout.write('\n'+'-'*N)
        sys.stdout.flush()
        all_numerical_x = []
        if True:
            slopes = d.get('slopes')
            all_numerical_x = v_th - slopes*dt
        else:
            mem_mon.insert_spikes(spike_mon, v_th)
            for spikes, V in zip(spike_mon.spiketimes.values(), mem_mon.values):
                figure(1)
                clf()
                nspikes = len(spikes)
                ISI = diff(spikes)
                mslope, slopes = nt.firing_slope(V, spikes)
                numerical_x = v_th - slopes*dt
                all_numerical_x = append(all_numerical_x, numerical_x)
        h, b = histogram(all_numerical_x, density=True, bins=50)
        x = linspace(min(all_numerical_x), max(all_numerical_x), 1000)
        dx = diff(x)[0]
        #t = duration
        #Mt = mu*tau*(1-exp(-t/tau))
        t = 0.1*ms
        Mt = mu*tau + (v_th - mu*tau)*exp(-t/tau)
        Vt = (((sigma**2) * tau)/2) * (1-exp((-2*t)/tau))
        fxt = ((2*pi*Vt)**-0.5) * exp(-((x-Mt)**2) / (2*Vt))
        h *= trapz(fxt, dx=dx)

        mu_str = display_in_unit(mu, volt/second)
        si_str = display_in_unit(sigma, volt/second**0.5)
        Mt_str = display_in_unit(Mt, volt)
        Vt_str = display_in_unit(Vt, volt**2)

        subplot2grid((2,2), (0,0))
        plot(b[:-1], h, color='blue')
        title("Numerical")

        subplot2grid((2,2), (0,1))
        plot(x, fxt, color='green')
        title("Theoretical")

        subplot2grid((2,2), (1,0), colspan=2)
        plot(b[:-1], h, label="Numerical", color='blue')
        plot(x, fxt, label="Theoretical", color='green')
        legend(loc='best')

        suptitle("Transition density mu=%s, sigma=%s" % (mu_str, si_str))
        savefig("%s_%i_trandens.png" % (filehead, progress))
        continue
        if nspikes > 20:
            '''First passage time density (ISI dist)'''
            ISIdist_h, ISIdist_b = histogram(ISI, density=True,
                    bins=30)
            t = linspace(0, ISIdist_b[-1], 1000)
            if regime == -1: # sub threshold
                lam = sigma*sqrt(pi*tau) / (v_th-mu*tau) *\
                        exp(-(v_th-mu*tau)**2 / (sigma**2 * tau))
                fptd = lam * exp(-lam*t)
                theta = (v_th-mu*tau) / (sigma*sqrt(tau))
                theta_actual = append(theta_actual, theta)

                mISI = mean(ISI)
                theta_rng = linspace(-10, 10, 10000)
                mISI_est = theta_rng/sqrt(pi) * exp(theta_rng**2)
                err = abs(mISI - mISI_est)
                best_theta = theta_rng[where(err == min(err))]
                theta_est = append(theta_est, best_theta)

            elif regime == 0: # near threshold
                exp2tt = exp(2*t/tau)
                fptd = 2*v_th*exp2tt / \
                        (sqrt(pi*tau**3*sigma**2) * (exp2tt-1)**(3./2)) * \
                        exp(-v_th**2 / (sigma**2 * tau * (exp2tt-1)))
                mu_est = append(mu_est, v_th/tau)
                mu_actual = append(mu_actual, mu)

                sigma_est = append(sigma_est, sqrt(
                        mean(
                            2*v_th**2 / (tau*(exp(2*diff(spikes)/tau)-1))
                            )
                        ))
                sigma_actual = append(sigma_actual, sigma)
            elif regime == 1: # super threshold
                exp2tt = exp(2*t/tau)
                fptd = 2*v_th*exp2tt / \
                        (sqrt(pi*tau**3*sigma**2) * (exp2tt-1)**(3./2)) * \
                        exp(-v_th**2 / (sigma**2 * tau * (exp2tt-1)))
                Z1 = mean(exp(ISI/tau))
                Z2 = mean(exp(2*ISI/tau))
                mu_est = append(mu_est,
                        v_th*Z1 / (tau*(Z1-1))
                        )
                mu_actual = append(mu_actual, mu)
                sigma_est = append(sigma_est, sqrt(
                    2*v_th**2 * (Z2 - Z1**2) / (tau * (Z2 - 1)*(Z1 - 1)**2)
                    ))
                sigma_actual = append(sigma_actual, sigma)
            else:
                print "Something went wrong: regime value is %i" % regime

            figure(1)
            clf()
            subplot2grid((2,2), (0,0))
            plot(ISIdist_b[:-1], ISIdist_h, color='blue')
            title("Numerical")
            subplot2grid((2,2), (0,1))
            plot(t, fptd, color='green')
            title("Theoretical")
            subplot2grid((2,2), (1,0), colspan=2)
            plot(ISIdist_b[:-1], ISIdist_h, label='Numerical', color='blue')
            plot(t, fptd, label="Theoretical", color='green')
            legend(loc='best')
            suptitle("First passage time density: %s\
                    \nmu=%s, sigma=%s, Nspikes=%i" % (regime_str,
                        mu_str, si_str, nspikes))
            subplots_adjust(top=0.8)
            savefig("%s_%i_fptd.png" % (filehead, progress))
        progress += 1
        sys.stdout.write('\r'+'|'*progress+'-'*(N-progress))
        sys.stdout.flush()
        i+=1

    mu_act_unique = sorted(unique(mu_actual))
    mu_est_averages = [mean(mu_est[where(m == mu_actual)[0]])
            for m in mu_act_unique]
    sigma_act_unique = sorted(unique(sigma_actual))
    sigma_est_averages = [mean(sigma_est[where(s == sigma_actual)[0]])
            for s in sigma_act_unique]
    theta_act_unique = sorted(unique(theta_actual))
    theta_est_averages = [mean(theta_est[where(th == theta_actual)[0]])
            for th in theta_act_unique]
    figure(2)
    if len(mu_est) > 0:
        clf()
        plot(mu_actual, mu_est, 'o', label="ind sims")
        plot(mu_act_unique, mu_est_averages, '-', color='blue',
                label="mean")
        line = linspace(float(min(mu_actual)), float(max(mu_actual)))
        plot(line, line, color='red', label="1:1")
        xlabel("actual")
        ylabel("estimated")
        title("Mu")
        legend(loc="best")
        savefig("tmp/mu_est.png")

    if len(sigma_est) > 0:
        clf()
        plot(sigma_actual, sigma_est, 'o', label="ind sims")
        plot(sigma_act_unique, sigma_est_averages, '-', color='blue',
                label="mean")
        line = linspace(float(min(sigma_actual)), float(max(sigma_actual)))
        plot(line, line, color='red', label="1:1")
        xlabel("actual")
        ylabel("estimated")
        title("Sigma")
        legend(loc="best")
        savefig("tmp/sigma_est.png")

    if len(theta_est) > 0:
        clf()
        plot(theta_actual, theta_est, 'o', label="ind sims")
        plot(theta_act_unique, theta_est_averages, '-', color='blue',
                label="mean")
        line = linspace(float(min(theta_actual)), float(max(theta_actual)))
        plot(line, line, color='red', label="1:1")
        xlabel("actual")
        ylabel("estimated")
        title("Theta")
        legend(loc="best")
        savefig("tmp/theta_est.png")


