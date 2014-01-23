from brian import *
from brian.library.random_processes import *
import neurotools
from sys import stdout
from multiprocessing import Process, Lock, Queue
import time
import pickle
from numpy import arange, random
import os

def ousim(in_mu,in_sigma,res,l):
    '''
    Single neuron simulation for parallel execution.
    One thread per neuron.
    '''
    duration = 1000*ms
    v_reset = 0*mV
    v_rest = 0*mV
    v_th = 15*mV
    tau_mem = 10*ms
    capa = 1*mfarad
    refr = 2*ms
    w = 2*ms
    if os.name=='posix':
        random.seed(int((time.time()+in_mu/amp*1000+in_sigma/amp*1000000)))
    eqs = Equations('dv/dt = (-v-v_rest)/tau_mem + inp/capa : volt')
    eqs+=OrnsteinUhlenbeck('inp',mu=in_mu,sigma=in_sigma,tau=tau_mem)
    nrn = NeuronGroup(1,eqs,threshold=v_th,reset=v_reset,refractory=refr)
    stdout.flush()
    nrn.rest()

    v_mon = StateMonitor(nrn,'v',record=True)
    out_mon = SpikeMonitor(nrn)
    run(duration)
    if (out_mon[0].size > 1):
        #(sta_avg, sta_std, sta_wins) = neurotools.sta(v_mon[0],out_mon[0],w,0.1*ms)
        #slope = mean(diff(sta_avg))
        simslope = diff(v_mon[0])
        pos_slope = simslope[simslope > 0]
        slope = mean(pos_slope)
        '''
        also try mean POSITIVE slopes
        and slope distribution
        '''

    else:
        slope = 0

    try:
        res.put([in_mu,in_sigma,slope])
    except Exception as inst:
        l.acquire()
        print "Exception of type",type(inst),"thrown"
        print "mu:",in_mu,"  sigma:",in_sigma,"  slope:",slope
        print " "
        l.release()
    res.close()
    sys.exit()



if __name__=='__main__':
    savelocation = "/home/achilleas/Documents/Uni/working_dir/simout.py/exp.tmp/"

    '''
    Ranges for mu and sigma, the input parameters.
    We are testing how the slope reacts to these two parameters
    '''

    sigma_0 = 0.*mA
    sigma_end = 15.*mA
    sigma_step = 5.*mA
    mu_0 = 0.*mA
    mu_end = 15.*mA
    mu_step = 5.*mA
    n_samples = 9

    n_threads = 8

    sigma_range = arange(sigma_0,sigma_end+(sigma_step/10),sigma_step)
    mu_range = arange(mu_0,mu_end+(mu_step/10),mu_step)
    numsims = len(mu_range)*len(sigma_range)*n_samples
    print "Preparing to run",numsims,"simulations comprised of\n",len(mu_range),"mu values\n",len(sigma_range),"sigma values and\n",n_samples,"samples of each."

    if (numsims > 300):
        print "\nWARNING! Running this many simulations (",numsims,") may cause loss of data or the program to hang due to problems with the multiprocessing Queue. It is suggested that you keep the number of simulations below 300."
        answer = ""
        while not answer in ("y", "Y", "n", "N"):
            answer = raw_input("Are you sure you want to continue [y/n]? ")

        if (answer == "n" or answer == "N"):
            print "Exiting."
            sys.exit(1)

    p_list = []
    running = []
    '''
    The Queue usually hangs between 260-300 entries. I think this is a limitation of the object itself. Perhaps we could populate the array directly or use multiple Queues
    '''
    results = Queue(numsims)
    lock = Lock()
    for mu in mu_range:
        for sigma in sigma_range:
            for n in arange(n_samples):
                p_list.append(Process(target=ousim, args=(mu*amp,sigma*amp,results,lock)))

    stdout.write("\n")
    nalive = numsims
    p_i = 0
    while (nalive > 0):
        nalive = 0
        for p in p_list:
            if p.is_alive():
                nalive+=1

        while (nalive < n_threads and p_i < numsims):
            p_list[p_i].start()
            p_i+=1
            nalive+=1

        stdout.write("{0:d} processes running. Total: {1:d}. Results: {2:d}/{3:d}    \r".format(nalive,len(p_list),results.qsize(),numsims))
        stdout.flush()

    stdout.write("\n")
    resarray = zeros([len(mu_range),len(sigma_range)])

    if (results.qsize() == numsims):
        print "Done! All results added to queue"
    else:
        print "WARNING! NOT ALL PROCESSES HAVE PROVIDED DATA"
        print numsims-results.qsize(),"results missing!!!"
        answer = ""
        while not answer in ("y", "Y", "n", "N"):
            answer = raw_input("Plot results anyway [y/n]? ")

        if answer in ("n", "N"):
            print "Exiting."
            sys.exit(1)


    while (not results.empty()):
        r = results.get()
        r_i = mu_range == r[0]
        r_j = sigma_range == r[1]
        resarray[r_i,r_j] += r[2]

    resarray = resarray/n_samples
    print resarray
    clf()
    contourf(sigma_range,mu_range,resarray,numsims)
    colorbar()
    xlabel("Sigma")
    ylabel("Mu")
    savefig(savelocation+sys.argv[0]+str(int(time.time())))
    show()
    stdout.write("Program run complete!\n")
    

