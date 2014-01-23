nspikes = []
duration = 1*second
xdata = multiply(add(sigma_offs, sigma_peaks), 1000)
xlabelstr = "$\sigma_o+\sigma_A$ (mV)"
for st in sts:
    nspikes.append(len(st[0])/duration)
figure()
subplot(2,1,1)
plot(xdata, multiply(mslopes, 1000), '.b', markersize=15)
ylabel("Slope (mV/s)")
subplot(2,1,2)
plot(xdata, nspikes, '.b', markersize=15)
ylabel("$f_{out}$ (Hz)")
xlabel(xlabelstr)
suptitle("Slope and $f_{out}$ Vs $\sigma$")
show()

figure()
xdata = multiply(add(mu_offs, mu_peaks), 1000)
xlabelstr = "$\mu_o+\mu_A$ (mV)"
subplot(2,1,1)
plot(xdata, multiply(mslopes, 1000), '.b', markersize=15)
ylabel("Slope (mV/s)")
subplot(2,1,2)
plot(xdata, nspikes, '.b', markersize=15)
ylabel("$f_{out}$ (Hz)")
xlabel(xlabelstr)
suptitle("Slope and $f_{out}$ Vs $\mu$")
show()

figure()
xdata = freqs
xlabelstr = "$f$ (Hz)"
subplot(2,1,1)
plot(xdata, multiply(mslopes, 1000), '.b', markersize=15)
ylabel("Slope (mV/s)")
subplot(2,1,2)
plot(xdata, nspikes, '.b', markersize=15)
ylabel("$f_{out}$ (Hz)")
xlabel(xlabelstr)
suptitle("Slope and $f_{out}$ Vs $f$ (sine wave frequency)")
show()
