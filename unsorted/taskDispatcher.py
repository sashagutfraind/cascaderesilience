'''
Provides an architecture for creating threads and running pp jobs


    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


'''

import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
import numpy as np
import numpy.random as npr
import random
import matplotlib
matplotlib.use('PS')
import matplotlib.pylab as pylab
import networkx as nx #must be after pylab to prevent problems
import time, os, sys
import pp
import threading, thread


defaultServerList = (
                       'localhost:60001',
                       'localhost:60002',
                       'localhost:60003',
                       'localhost:60004',
     #                 'localhost:60005',
     #                 'localhost:60006',
     #                 'localhost:60007',
     #                 'localhost:60008',
     #                 'localhost:60009',
     #                 'localhost:60010',

                     )
#print 
#print 'Warning: pp is disabled - set server list correctly'
#defaultServerList = ('localhost', )
#defaultServerList = tuple()
defaultSecret     = 'cheb'

#to start, use:
#cd <path to netstorm code>
#python <path to python scripts>/ppserver.py -r -s cheb [-w 4]

#useful information about problems and good coding patterns:
#http://effbot.org/zone/thread-synchronization.htm

class DictionaryStorageThread ( threading.Thread ):
   def __init__( self, fn, args, storage, result_key, storage_lock, debug=False):
      if not debug:
         threading.Thread.__init__ (self)
      self.fn          = fn
      self.args        = args
      self.storage     = storage
      self.result_key  = result_key
      self.storage_lock= storage_lock 
      self.debug       = debug

   def isAlive(self):
      if not self.debug:
         return threading.Thread.isAlive(self)
      else:
         return False

   def run(self):
      try:
          #print 'This is thread ' + str ( self ) + ' speaking.'
          ret = self.fn(self.args)
          self.storage_lock.acquire()
          self.storage[self.result_key] = ret
      except Exception, inst:
          print 'Error executing thread: '
          print inst
          raise
      finally:
          if self.storage_lock.locked:
              self.storage_lock.release()

   def start(self):
      if not self.debug:
         threading.Thread.start(self)
      else:
         self.run()


def create_job_server(simulatedPP=True, serverList=defaultServerList, secret=defaultSecret):
    #create a true or simulated job server

    #a useful trick for anoymous classes, used when no job servers are available:
    #http://norvig.com/python-iaq.html
    class Struct:
        def __init__(self, **entries): self.__dict__.update(entries)

    if simulatedPP:
        job_server = Struct()
        job_server.submit = lambda func,args=tuple(),depfuncs=tuple(),modules=tuple(),globals=tuple():lambda :func(*args) #creates a lambda fn which takes no arguments and returns the evaluated number  
        job_server.destroy = lambda : ''
        job_server.get_ncpus = lambda : -10
        job_server.get_active_nodes = lambda : {}
        job_server.print_stats = lambda : ''
        job_server.secret = 'epo20pdosl;dksldkmm'
        job_server.ppservers = [('simulated', 60000)]
        job_server.simulated   = True
        job_server.wait = lambda : ''
        print 'Created a SIMULATED job server.'
        print
    else:
        print 'Server list: \n' + str(serverList)
        job_server = pp.Server(ppservers=tuple(serverList), loglevel=pp.logging.DEBUG, restart=True, secret=secret)
        
        job_server.set_ncpus(0)  #making all jobs remote removes unexpected errors
        time.sleep(1)
        print 'Active nodes: \n' + str(job_server.get_active_nodes())
        
    return job_server


#used by testThreadBasic():
theTestGlobalVar = 1

def callMethod(inputFilename):
    import cPickle, os, sys
    import findBest
    
    print 'Loading '+inputFilename
    
    f = open(inputFilename, 'r') 
    executionParams = cPickle.load(f)
    f.close()
    
    fnName   = executionParams['fnFullName']
    paramSet = executionParams['paramSet']
    outputFilename = executionParams['outputFilename']
    
    print 'Starting: '+fnName
    try:
        fnCode = eval(fnName)
        ret    = fnCode(paramSet)

        f = open(outputFilename, 'w')
        cPickle.dump(ret, f)
        f.close()
        fDone = open(outputFilename+'.Done', 'w')
        fDone.close()
    except Exception, inst:
        print 'Something went wrong...'
        import traceback
        traceback.print_exc(file=sys.stdout)
        print inst
        raise 
        fDone = open(outputFilename+'.Failed', 'w')
        fDone.close()

def testThreadBasic():
#a very primitive example
    class MyThread ( threading.Thread ):
       def run ( self ):
          global theTestGlobalVar
          print 'This is thread ' + str ( theTestGlobalVar ) + ' speaking.'
          print 'Hello and good bye.'
          theTestGlobalVar = theTestGlobalVar + 1

    for x in xrange ( 20 ):
       MyThread().start()


def testFnCall():
#a pretty generic example of dispersing a problem to threads, with safety
    class MyThread2 ( threading.Thread ):
       def __init__( self, fn, args, storage, storage_key, storage_lock):
          threading.Thread.__init__ (self)
          self.fn          = fn
          self.args        = args
          self.storage     = storage
          self.storage_key = storage_key
          self.storage_lock= storage_lock

       def run(self):
          #print 'This is thread ' + str ( self ) + ' speaking.'
          try:
              ret = self.fn(self.args)
              self.storage_lock.acquire()
              print 'locked'
              self.storage[self.storage_key] = ret
          except Exception, inst:
              print 'Error executing thread: '
              print inst
              raise
          finally:
              if self.storage_lock.locked:
                  self.storage_lock.release()
                  print 'unlocked'

    def factorial(args):
        g = args['g']
        return reduce(lambda x,y: x*y, xrange(1, g+1), 1)

    goals = []
    for i in range(20):
        goals.append(npr.randint(0, 50))

    storage_lock = thread.allocate_lock()
    results = {}
    threads = {}
    for g in goals:
       th = MyThread2(fn=factorial, args={'g':g}, storage=results, storage_key=g, storage_lock=storage_lock)
       threads[g] = th
       #all have the lock now
    for th_key in threads:
       threads[th_key].start()

    while threads != {}:
        g, th = threads.items()[0]
        if th.isAlive():
           print g, th
           print 'waiting...'
           time.sleep(1)
        else:
           print g, th
           threads.pop(g)

    print 'Printing result...'
    print results

def testThreads():
    testThreadBasic()
    testFnCall()
    
if __name__ == '__main__':
    import sys, getopt, re
    opts, args = getopt.getopt(sys.argv[1:], 'tP', ['test', 'Pkl='])

    inputFilename= None

    for o, a in opts:
       if o in ('-P', '--Pkl'):
          callMethod(inputFilename=a)
       if o in ('-t', '--test'):
          testThreads()


