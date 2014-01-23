import numpy as np
import sys
sys.path.append('../../QuickGA/')
import QuickGA as QGA

filename = 'mclownoise.npz'
archive = np.load(filename)
mu_offs_actual = archive['mu_offs_actual']
mu_t_est = archive['mu_t_est']

def fitnessfunc(individual, args):
    actual, mu_t_est = args
    weights = individual.chromosome
    est = np.array([np.average(mt, weights=weights) for mt in mu_t_est])
    error = abs(est-actual)/actual
    individual.fitness = np.mean(error)

optimiser = QGA.GA(20, 10, outputfile='weightopt.txt')
optimiser.fitnessfunc = fitnessfunc
goodind = optimiser.Individual(np.array([2.725634, -0.344084, 0.852556, 0.218468, 0.125677,
    0.098427, 1.455072, 0.333759, -0.425731, 0.858043]), 0)
optimiser.population[0] = goodind
optimiser.optimise(1000000, mu_offs_actual, mu_t_est)


