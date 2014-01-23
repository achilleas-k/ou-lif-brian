import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
from brian import *

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True) # use latex

'''
Plot freq est err vs mu peak est err
'''
filename = 'lotsofdata.npz'
archive = np.load(filename)
mu_peaks_actual = archive['mu_peaks_actual']
mu_peaks_est_bad = archive['mu_peaks_est_bad']
nzind = flatnonzero(mu_peaks_actual)
mu_peaks_actual = mu_peaks_actual[nzind]
mu_peaks_est_bad = array([mu_peaks_est_bad[n] for n in nzind])
fest_err = linspace(-0.1, 0.1, 21)
mu_peaks_est_err = array([mean(abs(mu_peaks_actual-mpeb)/mu_peaks_actual) for
                                mpeb in mu_peaks_est_bad.transpose()])
figure(num=1, figsize=(8, 6), dpi=100)
clf()
plot(fest_err, mu_peaks_est_err, 'k')
xlabel(r'$ \varepsilon_f $', size=15)
ylabel(r'$ \varepsilon_{\mu_p} $', size=15)
xticks(linspace(-0.1, 0.1, 11), size=15)
yticks(size=15)
axis(xmin=-0.1, xmax=0.1)
savefig('brain_research_figures/mu_err_vs_freq_err.png')
savefig('brain_research_figures/mu_err_vs_freq_err.pdf')
print "Figure saved!"

