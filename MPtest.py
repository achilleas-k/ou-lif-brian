from multiprocessing import Queue, Process
import os
from time import time, sleep

def afunc(snum):
    #sleep(3) # arbitrary run time
    aqueue.put(snum)
    return


plist = [] #remaining process list
rlist = [] #running process list
aqueue = Queue(10000)
maxprocs = 10
running = 0
finished = 0
remaining = 0

for snum in range(10000):
    newproc = Process(target=afunc, args=(snum,))
    plist.append(newproc)
    remaining += 1

while plist or rlist:
    while len(rlist) < maxprocs and plist:
        startproc = plist.pop()
        startproc.start()
        rlist.append((startproc, time()))
        #print "Process %i started" % (startproc.pid)
        remaining -= 1
        running += 1

    for rproc, stime in rlist[:]:
        if not rproc.is_alive():
            rlist.remove((rproc, stime))
            #print "Process %i finished after %f secs" % (
            #        rproc.pid, time()-stime)
            running -= 1
            finished += 1
        elif time()-stime > 60:
            print("Process %i has been running for %f secs" % (
                    rproc.pid, time()-stime))
    print("%i running, %i finished, %i remaining" % (
            running, finished, remaining))
    sleep(1)

aqueue.put('done')
for q in iter(aqueue.get, 'done'):
    print(q)

