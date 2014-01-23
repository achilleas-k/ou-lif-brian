def multi_scatter(data_items):
    h = 0
    for di in data_items:
        st = di[6][0]
        so = di[2]
        sa = di[3]
        scatter(st, ones(len(st))*h, label='h=%i, \sigma_o=%s, \sigma_A=%s' %
                (h, so*volt, sa*volt))
        h+=1
    if h > 0:
        axis([0, 2, -0.5, h-0.5])
    yticks([])

def multi_hist(data_items):
    combined_st = array([])
    for di in data_items:
        st = di[6][0]
        combined_st = append(combined_st, st)
    hist(combined_st, linspace(0, 2, 101))


# separate results into three categories based on noise type (none,
# constant, sinusoidal)
# ---
# also pick data which correspond to specific mu pair
# ---

des_mo = 0.010; des_ma = 0.015
des_s = 0.002
const_sigma = []
sin_sigma = []
no_sigma = []
no_sigma_items = []
const_sigma_items = []
sin_sigma_items = []
for d in data.itervalues():
    mo, ma, so, sa = d[0:4]
    freq = d[4]
    if so == 0 and sa == 0:
#        no_sigma.append(d)
        if mo == des_mo and ma == des_ma:
            no_sigma_items.append(d)
    elif so > 0 and sa == 0:
#        const_sigma.append(d)
        if mo == des_mo and ma == des_ma:
            const_sigma_items.append(d)
    elif so > 0 and sa > 0:
#        sin_sigma.append(d)
        if mo == des_mo and ma == des_ma:
            sin_sigma_items.append(d)

freq_ang = freq*Hz*2*pi
rc('text', usetex=True)


raster_fig = figure()
subplot(411)
title("No noise")
h = 0
multi_scatter(no_sigma_items)
subplot(412)
title("Constant noise")
multi_scatter(const_sigma_items)
subplot(413)
title("Sinusoidal noise")
multi_scatter(sin_sigma_items)
subplot(414)
title("Signal")
t = arange(0, 2.01, 0.01)
plot(t,sin(t*freq_ang))
yticks([])
axis([0,2,-1,1])
xlabel("t (sec)")

# Instead of plotting spike rasters, consider plotting the sum of the
# spike densities. This should show a curve that resembles the sine wave
# and becomes increasingly smooth for the three cases.

# This could either be a histogram with an appropriate bin-width, or the
# result of summing the spike trains and convolving with a gaussian
# kernel.
hist_fig = figure()
subplot(411)
title("No noise")
h = 0
multi_hist(no_sigma_items)
subplot(412)
title("Constant noise")
multi_hist(const_sigma_items)
subplot(413)
title("Sinusoidal noise")
multi_hist(sin_sigma_items)
subplot(414)
title("Signal")
t = arange(0, 2.01, 0.01)
plot(t,sin(t*freq_ang))
yticks([])
axis([0,2,-1,1])
xlabel("t (sec)")


# Contour (3D) plot with:
# X - time
# Y - mu_offs (or mu_amp; keep the other one fixed)
# Z - spike density (colour)
# This will show how the spike density increases with the signal
# amplitude. If we plot the same figure for constant and sinusoidal
# noise, we would see if the signal can be determined for a wider range
# of signal amplitudes in the case of sinusoidal noise. This should
# start from sub-threshold (very low amplitudes) and finish at
# super-threshold inputs to cover the entire range. The no-noise
# spiking will be included as reference.
# sigma_offs and sigma_amp should be fixed for all simulations
# additionally, in the case of constant noise, we should use
# sigma_offs+sigma_amp, as that is the maximal value that the sinusoidal
# noise has.

