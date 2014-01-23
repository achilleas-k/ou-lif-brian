ncols = 4
nrows = len(sigma_peaks)/ncols
plot_num = 0
subplots_adjust(hspace=1, wspace=0.5)
sort_index = sorted(range(len(sigma_peaks)), key=sigma_peaks.__getitem__)
sp_sorted = array(sigma_peaks)[sort_index]
slopes_sorted = array(slopes)[sort_index]
for sigm, slp in zip(sp_sorted, slopes_sorted):
    if sigm == 0:
        continue
    plot_num += 1
    subplot(nrows, ncols, plot_num)
    s = [si*1000 for si in flatten(slp)]
    srange = max(s) - min(s)
    binwidth = srange/10
    y, binedges = np.histogram(s, normed=True)
    bincentres = 0.5*(binedges[1:]+binedges[:-1])
    bar(binedges[:-1], y/sum(y), width=binwidth)
    plot(bincentres, y/sum(y), 'r', linewidth=3)
    xlabel('Slope (mV/s)')
    title('$\sigma_A$: %2.f mV (mean: %.1f)' %
            ((sigm*1000), round(mean(s), 1)))
    ax = axis()
    axis([ax[0],ax[1],0,0.4])
    yticks(arange(0, 0.41, 0.1))
    xticks(arange(0, ax[1]+0.001, ax[1]/5).round())
#    tick_params(direction='out')
#    savefig("slopedist_sigma_%s.png" % sigm)
#    clf()
suptitle('Slope distribution for each $\sigma_A$ value')
