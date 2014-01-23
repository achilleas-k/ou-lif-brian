pltnum = 0
allbins = zeros(10)
for ind in range(len(mslopes))[-20:]:
    f = freqs[ind]
    st = sts[ind]
    mu_peak = mu_peaks[ind]
    sigma_peak = sigma_peaks[ind]
    mu_off = mu_offs[ind]
    sigma_off = sigma_offs[ind]
#    pltnum+=1
#    subplot(4, 5, pltnum)
    left, bins = nt.spike_period_hist(st.values(), f, 1*second)
    allbins = allbins + bins

bar(arange(0,1,0.1), allbins, width=0.1)
axis([0,1,0,max(allbins)])
t = arange(0,1,1./100)
plot(t, max(allbins)/2+max(allbins)/2*sin(t*2*pi), 'r')
#    xticks([0, 1])
#    yticks([0, max(bins)])
#    title("$\mu_A$: %.2f, $\mu_o$: %.2f,\n $\sigma_A$: %.2f, $\sigma_o$: %.2f,\
#f: %i" %
#            (mu_peak, mu_off, sigma_peak, sigma_off, f),
#            {'fontsize': 'small'})

#subplots_adjust(top=0.96, bottom=0.02, left=0.01, right=0.99,
#            wspace=0.4, hspace=2)
show()
