import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
from brian import *

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True) # use latex

'''
Plot CV vs freq est error and histogram
'''
filename = 'lotsofdata.npz'
archive = np.load(filename)
CV = archive['CV']
freq_est = archive['freq_est']
freq_actual = archive['freq_actual']
ferr = abs(freq_est-freq_actual)/freq_actual


figure(num=1, figsize=(8, 12), dpi=100)
clf()

subplot(2,1,1)
scatter(CV, ferr, c='k')
ylabel(r'$ \varepsilon_f $', size=15)
yticks(size=15)
axis([0, 2.5, -0.1, 10])

errcv = CV[ferr > 0]
subplot(2,1,2)
hist(errcv, bins=10, normed=True, color='gray')
xlabel(r'CV', size=15)
axis([0, 2.5, 0, 10])

savefig('brain_research_figures/CV_freq_est_err.png')
savefig('brain_research_figures/CV_freq_est_err.pdf')
print "Figure saved!"

