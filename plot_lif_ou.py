from brian import *
import matplotlib as mpl

duration = 1000*ms
tau = 10*ms
V_th = 10*mV
t_refr = 0*ms
V_reset = 0*mV
V0 = 0*mV

mu_amp = 1*mV/ms
mu_offs = 1*mV/ms
sigma_amp = 0.5*mV/sqrt(ms)
sigma_offs = 0.5*mV/sqrt(ms)
freq = 20*Hz
freq_ang = freq*2*pi

dt = defaultclock.dt

def ousim():
    clock.reinit_default_clock()
    eqs =Equations('dV/dt = mu-(V+V0)/tau + sigma*I/sqrt(dt) : volt')
    eqs+=Equations('dI/dt = -I/dt + xi/sqrt(dt) : 1')
    eqs+=Equations('mu = mu_amp*sin(t*freq_ang) + mu_offs : volt/second')
    eqs+=Equations('sigma = sigma_amp*sin(t*freq_ang) + sigma_offs :'
                                                        ' volt/sqrt(second)')
    eqs.prepare()
    ounrn = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr,
                                                                reset=V_reset)
    ounrn.V = V0
    mu_mon = StateMonitor(ounrn, 'mu', record=True)
    sigma_mon = StateMonitor(ounrn, 'sigma', record=True)
    I_mon = StateMonitor(ounrn, 'I', record=True)
    V_mon = StateMonitor(ounrn, 'V', record=True)
    st_mon = SpikeMonitor(ounrn)
    eqs.prepare()
    run(duration)
    V_mon.insert_spikes(st_mon, value=V_th*2)
    times = V_mon.times
    membrane = V_mon[0]
    input_trace = (I_mon[0]*sigma_mon[0]/sqrt(dt)+mu_mon[0])
    return times-times[-1]*0.2, membrane, input_trace

def lifsim():
    clock.reinit_default_clock()
    eqs = Equations('dV/dt = (-V+V0)/tau : volt')
    lifnrn = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr,
                                                                reset=V_reset)
    pulse_times = arange(1.25, 20, 1)/freq
    pulse_spikes = []
    for pt in pulse_times:
        try:
            pp = PulsePacket(t=pt*second, n=100, sigma=1/(freq*8))
            pulse_spikes.extend(pp.spiketimes)
        except ValueError:
            print "Skipping pulse packet at %s" % (pt*second)
    pulse_input = SpikeGeneratorGroup(100, pulse_spikes)
    pulse_conn = Connection(pulse_input, lifnrn, 'V', weight=0.5*mV)
    poiss_input = PoissonGroup(100, freq)
    poiss_conn = Connection(poiss_input, lifnrn, 'V', weight=0.0*mV)
    V_mon = StateMonitor(lifnrn, 'V', record=True)
    pulse_mon = PopulationRateMonitor(pulse_input, bin=0.1*ms)
    poiss_mon = PopulationRateMonitor(poiss_input, bin=0.1*ms)
    st_mon = SpikeMonitor(lifnrn)
    eqs.prepare()
    run(duration)
    V_mon.insert_spikes(st_mon, value=V_th*2)
    times = V_mon.times
    membrane = V_mon[0]
    input_trace = (pulse_mon.rate + poiss_mon.rate)
    input_trace = movavg(input_trace, 50)
    input_trace = append(input_trace, zeros(len(times)-len(input_trace)))/50
    return times-times[-1]*0.2, membrane, input_trace

if __name__=='__main__':
    T_ou, V_ou, I_ou = ousim()
    print "OU done"

    T_lif, V_lif, I_lif = lifsim()
    print "LIF done"

    figure(num=1, figsize=(13.66,7.68), dpi=100)
    ax_limits = [0*second, duration*0.5, 0*mvolt, V_th*2]

    frows = 2
    fcols = 1
    fnum = 1

    #subplot(frows, fcols, fnum)
    #plot(T_ou, I_ou)
    #xlabel('t (sec)')
    #ylabel('I (volt)')
    #title('OU Input')
    #axis(ax_limits)
    #fnum += 1

    subplot(frows, fcols, fnum)
    plot(T_ou, V_ou)
    xlabel('t (sec)')
    ylabel('x (volt)')
    title('OU Voltage')
    axis(ax_limits)
    fnum += 1

    #subplot(frows, fcols, fnum)
    #plot(T_lif, I_lif)
    #xlabel('t')
    #ylabel('I')
    #title('LIF Input')
    #axis(ax_limits)

    subplot(frows, fcols, fnum)
    plot(T_lif, V_lif)
    xlabel('t (sec)')
    ylabel('V (volt)')
    title('LIF Voltage')
    axis(ax_limits)
    fnum += 1

    subplots_adjust(hspace=0.4)
    mpl.rcParams['font.size'] = 12
    savefig('ou_vs_lif.eps')
