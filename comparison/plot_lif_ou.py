from __future__ import print_function
from brian import *
import matplotlib as mpl
from spikerlib.metrics import kreuz

duration = 1000*ms
tau = 10*ms
V_th = 10*mV
V_th = 100*mV
t_refr = 0*ms
V_reset = 0*mV
V0 = 0*mV

mu_amp = 1.0*mV/ms
mu_offs = 1.0*mV/ms
sigma_amp = 0.1*mV/sqrt(ms)
sigma_offs = 0.1*mV/sqrt(ms)
freq = 20*Hz
freq_ang = freq*2*pi

dt = defaultclock.dt

def ousim():
    print("Setting up OU LIF simulation...")
    ounet = Network()
    clock.reinit_default_clock()
    eqs =Equations('dV/dt = mu-(V+V0)/tau + sigma*I/sqrt(dt) : volt')
    eqs+=Equations('dI/dt = -I/dt + xi/sqrt(dt) : 1')
    eqs+=Equations('mu = mu_amp*sin(t*freq_ang) + mu_offs : volt/second')
    eqs+=Equations('sigma = sigma_amp*sin(t*freq_ang) + sigma_offs :'
                                                        ' volt/sqrt(second)')
    eqs.prepare()
    ounrn = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr,
                                                                reset=V_reset)
    ounet.add(ounrn)

    ounrn.V = V0
    print("mu_amp:    {} mV/ms,       mu_offs:    {} mV/ms\n"
          "sigma_amp: {} mV/sqrt(ms), sigma_offs: {} mV/sqrt(ms)\n"
          "frequency: {} Hz".format(
              mu_amp, mu_offs, sigma_amp*sqrt(ms)/mV, sigma_offs*sqrt(ms)/mV,
              freq))
    print("Configuring monitors...")
    mu_mon = StateMonitor(ounrn, 'mu', record=True)
    sigma_mon = StateMonitor(ounrn, 'sigma', record=True)
    I_mon = StateMonitor(ounrn, 'I', record=True)
    V_mon = StateMonitor(ounrn, 'V', record=True)
    st_mon = SpikeMonitor(ounrn)
    ounet.add(mu_mon, sigma_mon, I_mon, V_mon, st_mon)

    print("Running simulation for {} seconds".format(duration))
    ounet.run(duration)
    print("OU done")

    V_mon.insert_spikes(st_mon, value=V_th*2)
    times = V_mon.times
    membrane = V_mon[0]
    input_trace = (I_mon[0]*sigma_mon[0]/sqrt(dt)+mu_mon[0])
    return times-times[-1]*0.2, st_mon, membrane, input_trace

def lifsim():
    print("Setting up LIF simulation...")
    lifnet = Network()
    clock.reinit_default_clock()
    eqs = Equations('dV/dt = (-V+V0)/tau : volt')
    eqs.prepare()
    lifnrn = NeuronGroup(1, eqs, threshold=V_th, refractory=t_refr,
                                                                reset=V_reset)
    lifnet.add(lifnrn)
    pulse_times = arange(1.25, 20, 1)/freq
    pulse_spikes = []
    print("Generating input spike trains...")
    Nin = 2000
    sigma = 1/(freq*5)
    weight = 50*mV/Nin
    for pt in pulse_times:
        try:
            pp = PulsePacket(t=pt*second, n=Nin, sigma=sigma)
            pulse_spikes.extend(pp.spiketimes)
        except ValueError:
            print("Skipping pulse packet at %s" % (pt*second))
    pulse_input = SpikeGeneratorGroup(Nin, pulse_spikes)
    pulse_conn = Connection(pulse_input, lifnrn, 'V', weight=weight)
    # poiss_input = PoissonGroup(100, freq)
    # poiss_conn = Connection(poiss_input, lifnrn, 'V', weight=0.0*mV)
    lifnet.add(pulse_input, pulse_conn)

    print("N_in: {}, Npulses: {}\n"
          "S_in: {}, sigma_in: {} ms, weight: {} mV".format(
              Nin, len(pulse_times), 1.0, sigma*1000, weight*1000))

    print("Configuring monitors...")
    V_mon = StateMonitor(lifnrn, 'V', record=True)
    pulse_mon = PopulationRateMonitor(pulse_input, bin=0.1*ms)
    # poiss_mon = PopulationRateMonitor(poiss_input, bin=0.1*ms)
    st_mon = SpikeMonitor(lifnrn)
    lifnet.add(V_mon, pulse_mon, st_mon)#, poiss_mon, st_mon)

    print("Running simulation for {} seconds".format(duration))
    lifnet.run(duration)
    print("LIF done")

    V_mon.insert_spikes(st_mon, value=V_th*2)
    times = V_mon.times
    membrane = V_mon[0]
    input_trace = (pulse_mon.rate) # + poiss_mon.rate)
    input_trace = movavg(input_trace, 50)
    input_trace = append(input_trace, zeros(len(times)-len(input_trace)))/50
    return times-times[-1]*0.2, st_mon, membrane, input_trace

if __name__=='__main__':
    if (V_th > 20*mV):
        print("SUB-THRESHOLD SIMULATIONS")
        fnamesuffix = "nosp"
    else:
        print("SUPRA-THRESHOLD SIMULATIONS")
        fnamesuffix = "sp"
    T_ou, S_ou, V_ou, I_ou = ousim()
    print("\n\n")
    T_lif, S_lif, V_lif, I_lif = lifsim()

    ax_limits= [0, 500, 0, 20]

    figure(figsize=(8, 6))
    subplot2grid((5, 1), (0, 0), rowspan=4, colspan=1)
    plot(T_lif*1000, V_lif*1000)
    plot(T_ou*1000, V_ou*1000)
    axis(ax_limits)
    ylabel("mV")
    xt, _ = xticks()
    xticks(xt, [])
    subplot2grid((5, 1), (4, 0), rowspan=1, colspan=1)
    plot(T_ou*1000, abs(V_ou-V_lif)*1000)
    xlabel("t (ms)")
    ylabel("mV")
    yticks([0, 2.5, 5])
    axis(ax_limits)
    axis(ymax=5)
    mpl.rcParams["font.size"] = 12
    subplots_adjust(left=0.1, top=0.95, bottom=0.1, right=0.95, hspace=0.2)
    savefig("ou_vs_lif_"+fnamesuffix+".pdf")

    figure(figsize=(8, 3))
    plot(T_ou*1000, V_ou*1000)
    xlabel("t (ms)")
    ylabel("mV")
    axis(ax_limits)
    mpl.rcParams["font.size"] = 12
    subplots_adjust(left=0.1, top=0.95, bottom=0.2, right=0.95)
    savefig("ou_sin_"+fnamesuffix+".pdf")

    figure(figsize=(8, 3))
    plot(T_lif*1000, V_lif*1000)
    xlabel("t (ms)")
    ylabel("mV")
    axis(ax_limits)

    mpl.rcParams["font.size"] = 12
    subplots_adjust(left=0.1, top=0.95, bottom=0.2, right=0.95)
    savefig("lif_sin_"+fnamesuffix+".pdf")

    dist = kreuz.distance(S_lif.spiketimes[0], S_ou.spiketimes[0],
                          0*second, duration, duration/(1*ms))
    kdist = np.trapz(dist[1], dist[0])
    print("Spike train distance: {}".format(kdist))
    memdiff = np.trapz(V_lif, T_lif)-np.trapz(V_ou, T_ou)
    print("Mem potential diff  : {}".format(memdiff))
    # show()
