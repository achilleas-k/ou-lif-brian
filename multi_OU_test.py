from brian import *
from brian.library.random_processes import *
from multiprocessing import Process, Queue
import matplotlib.pyplot as plt
import time

def ousim(in_mu,in_sigma,spikes,mems):
    '''
    Single neuron simulation for parallel execution.
    One thread per neuron.
    '''
    
    '''
    Reseeding numpy rng for parallel execution
    '''
    np.random.seed(int(time.time()*1000))

    duration = 5000*ms
    v_reset = 0*mV
    v_rest = 0*mV
    v_th = 15*mV
    tau_mem = 10*ms
    capa = 1*mfarad
    refr = 2*ms
    eqs = Equations('dv/dt = (-v-v_rest)/tau_mem + inp/capa : volt')
    eqs+=OrnsteinUhlenbeck('inp',mu=in_mu,sigma=in_sigma,tau=tau_mem)
    nrn = NeuronGroup(1,eqs,threshold=v_th,reset=v_reset,refractory=refr)
    nrn.rest()
    out_mon = SpikeMonitor(nrn)
    mem_mon = StateMonitor(nrn,'v',record=True)
    print "Starting sim m:",in_mu,", s:",in_sigma
    run(duration)
    spikes.put(out_mon)
    mems.put(mem_mon)    
    print "Finished sim m:",in_mu,", s:",in_sigma
    sys.exit()

if __name__=='__main__':
    n_samples = 4    
    spikes = Queue()
    mems = Queue()
    p_list = []
    pid = 0
    mu = 0.001*amp
    sigma = 0.001*amp
    for n in arange(n_samples):        
        p_list.append(Process(target=ousim, args=(mu,sigma,spikes,mems)))
        p_list[-1].start()
   
    mem_results = []
    out_results = []
    while len(mem_results) < n_samples or len(out_results) < n_samples:
        mem_results.append(mems.get())
        out_results.append(spikes.get())
        print "res size [m,o]:",len(mem_results),len(out_results)
        time.sleep(0.1)        

    plt.hold(True)
    for i in arange(len(mem_results)):
        '''
        mon = spikes.get()
        print mon[0]
        '''
        print out_results[i].nspikes
        mem = mem_results[i]
        print "mem m:",mean(mem[0]),"s:",std(mem[0])
#        plt.plot(mem.times,mem[0])

#    plt.legend(("1","2","3","4"))
#    plt.show()
    sys.exit(0)

