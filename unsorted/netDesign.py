'''
Network Generation Program
    - base class, constant design

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


Note:
    We use "cost" to refer to fitness (historical reasons)
'''

import sys, os
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
sys.path.append('/nh/nest/u/gfriend/lib/python')
import networkx as nx
import numpy as np
import numpy.random as npr
import time, pdb
import random
import pp

import resMetrics
import effMetrics
import graphStats


def defaultSIRcost(G, costParamsInput, fullData=True):
    sys.path.append('/home/gfriend/lib/python')
    sys.path.append('/homes/gfriend/lib/python')
    curWorkDir = costParamsInput.get('curWorkDir', None)
    if curWorkDir != None:
        os.chdir(curWorkDir)
    costParams = costParamsInput.copy()
    if 'tau' not in costParamsInput:
        costParams['tau'] = 0.5
        print 'Warning: using detault tau (%f)!'%costParams['tau']
    nn = G.number_of_nodes()
    if not costParams.get('weighted', False):
        avgEpidemic, stdevEpidemic = resMetrics.extentSIRSimX(G, params=costParams)
    else:
        avgEpidemic, stdevEpidemic = resMetrics.extentSIRSim(G, params=costParams)
    efficiency = effMetrics.ringReach(G, costParams, normalize=True)

    #fixme: consider correcting for self

    r  = costParams.get('r', 0.5)
    if 'r' not in costParamsInput:
        print 'Warning: using detault r (%f)!'%r
    if nn > 1: 
        resilience = 1.0 - (avgEpidemic-1.0) / (nn-1.0)
        #the fraction of the network (excluding the starting node) that is not infected
    else:
        resilience = 1.0
    fitness    = r*resilience + (1.0-r)*efficiency

    #print 'defaultSIRcost completed computation!'
    if fullData:
        return {'resilience':float(resilience), 'efficiency':float(efficiency), 'fitness':float(fitness)}
    else:
        return float(fitness)


def multipleNetsCostHelper(nets, costFn, costParamsInput, fullData=True):
    costs = []
    for netIdx, G in enumerate(nets):
        c = costFn(G=G, costParamsInput=costParamsInput, fullData=fullData)
        costs.append(c)

    return costs

class baseDesign:
    nets      = [] #network1, network2 etc. #Graph([0])
    netReports= []
    cost      = np.inf  #current. this should never be set to None -> causes failure when comparing
    fitness   = None

    costParams   = {} #for the cost fn
    fixedParams  = {} #immutate, like nn
    configParams = {} #mutable

    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        if x != None:
            self.nets        = x.nets[:]
            self.netReports  = x.netReports[:] #todo: maybe too shallow?
            self.cost        = x.cost


            self.netReports = []
            for report in x.netReports:
                self.netReports.append(report.copy())

            costParams    = x.costParams
            fixedParams   = x.fixedParams
            configParams  = x.configParams
            
        if costParams != None:
            self.costParams = costParams.copy()
        else:
            self.costParams = {'costFn': defaultSIRcost, 'tau':0.5, 'r':0.5, 'g':1.0,
                               'minReplications': 40, 'tolerance': 0.01, }

        if fixedParams != None:
            self.fixedParams = fixedParams.copy()
        else:
            self.fixedParams = {'nn': 100, 'numNets': 5, 'outputLevel':set(['debug'])}

        if configParams != None:
            self.configParams = configParams.copy()
        else:
            self.configParams = {} 

        #self.costParams['costFn'] = globals().get(self.costParams['costFnName'], defaultSIRcost)
        #self.costParams['costFn'] = locals()[self.costParams['costFnName']]

    def __str__(self):
        return '\nnets=%s\ncost=%s\nconfigParams=%s\ncostParams=%s\nfixedParams=%s\n'%(str(self.nets), str(self.cost), str(self.configParams), str(self.costParams), str(self.fixedParams))

    def __repr__(self):
        return self.__str__()
        
    def buildNet(self):
        #override me!
        raise Exception, 'call of abstract method'
        #G = nx.Graph()
        #G.add_edge(0,1)
        #return G

    def buildNets(self, seeds=None, updateInternals=True):
        #creates networks: updates internals, but could also be used to generate instances from this configuration 
        if seeds == None:
            seeds = []
            for i in xrange(self.fixedParams['numNets']):
                #s = int(time.time()*100 % 10000)
                seeds.append( int(npr.randint(1E6)) )
        nets       = []
        netReports = []
        for s in seeds:
            npr.seed(s)
            random.seed(s)
            ###
            net = self.buildNet() #sent to subclass #todo: this could potentially be dispatched to job_server
            ###
            #if net.selfloop_edges() != []: 
            #    net.delete_edges_from(net.selfloop_edges())
            #    if self.outputState('debug'): print 'self-loop edges detected and removed in '+str(self.__class__)
                #raise RuntimeError, 'Networks have self-loop edges!!!'
            nets.append(net)
            netReports.append({'seed':s})
            #netReports.append({'seed':s, 
            #                   'avgDegree':graphStats.avgDegree(net), 
            #                   's_metric': graphStats.sMetric(net)})
            #those reports are a bad idea - most networks will be thrown away in any case.  
            #better compute them at the end, by recreating the networks
        if updateInternals: #todo: some of this data could be made private
            self.nets       = nets
            self.netReports = netReports
        return nets, netReports

    def cross(self, parent2):
        childConfig = {}
        for param in self.configParams:
            if npr.rand() > 0.5:
                childConfig[param] = self.configParams[param]
            else:
                childConfig[param] = parent2.configParams[param]
         
        return self.__class__(fixedParams=self.fixedParams, 
                              costParams=self.costParams, 
                              configParams=childConfig)


    def getParam(self, paramName):
        if paramName in self.fixedParams:
            return self.fixedParams[paramName]
        if paramName in self.costParams:
            return self.costParams[paramName]
        if paramName in self.configParams:
            return self.configParams[paramName]
        raise ValueError, 'Unknown parameter: %s'%paramName

    def evaluate(self, configParams=None, nets=None, seeds=None):
        #dispatches costing through pp.  to achieve any speedup must be combined with threading on e.g. tau scan
        self.purge()
        if configParams != None: self.configParams = configParams.copy()  #since it will update self.cost etc

        if self.costParams.get('analytic'):
            return self.evaluateAnalytic()
        
        if nets == None and seeds == None:
            self.buildNets()
        elif nets != None and seeds == None:
            self.nets       = nets
            self.netReports = [{} for net in nets]
        elif nets == None and seeds != None:
            self.buildNets(seeds=seeds)
        else:
            self.nets       = nets
            self.netReports = [{} for net in nets]
            for i,seed in enumerate(seeds):
                self.netReports[i]['seed']=seed

        job_server = self.costParams.get('job_server', None)
        if self.outputState('debug'): print 'Evaluating... job_server: ' + str(job_server)
        
        costs  = [] 
        costFn = self.costParams['costFn']
        if job_server != None:
            prunedParams = self.costParams.copy()
            prunedParams.pop('job_server')  #cannot be pickled
            prunedParams['curWorkDir'] = os.getcwd()
            ppParams     = (self.nets, costFn, prunedParams) 
            job          = job_server.submit(func=multipleNetsCostHelper, args=ppParams, depfuncs=tuple(), 
                                             modules=tuple(['networkx', 'netDesign', #NX must be first
                                                             'effMetrics','resMetrics', 'os' ]), globals=tuple())
            ret = job()
        else:
            ret = multipleNetsCostHelper(self.nets, costFn, self.costParams)
        
        if type(ret) != type(None):  #success
            if self.outputState('debug'):
                print 'Job returned successfully.'
                print 'ret=\n%s'%str(ret)
            for netNum,netCost in enumerate(ret):
                 self.netReports[netNum]['cost'] = netCost
        else:
            print 'Error!'
            print type(ret)
            print 'baseDesign.evaluate()'
            print 'configParams: \n' + str(self.configParams)
            print 'costParams:   \n' + str(self.costParams)
            print 'fixedParams:  \n' + str(self.fixedParams)
            raise RuntimeError, 'Error: cost evaluation job failed!'

        if len(self.netReports) == 0:
            print 'WARNING: no networks generated!'
            return 0.0
        
        self.fitnessStats()
        return self.cost


    def evaluate_old(self, configParams=None, nets=None):
        #obsolete version of evaluate: building and costing is done locally
        self.purge()
        if configParams != None: self.configParams = configParams.copy()  #since it will update self.cost etc
        job_server = self.costParams.get('job_server', None)
        if self.outputState('debug'): print 'Evaluating... job_server: ' + str(job_server)
        if nets == None:
            self.buildNets()
        else:
            self.nets       = nets
            self.netReports = None

        jobs   = []
        costs  = [] 
        costFn = self.costParams['costFn']
        for netIdx, G in enumerate(self.nets):
            if job_server == None:
                c   = costFn(G, self.costParams)
                job = self.evaluateHelper(c)
            else:
                prunedParams = self.costParams.copy()
                prunedParams.pop('job_server')  #cannot be pickled
                prunedParams['curWorkDir'] = os.getcwd()
                ppParams     = (G, prunedParams) 
                job          = job_server.submit(func=costFn, args=ppParams, depfuncs=tuple(), 
                                                 modules=tuple(['networkx', 'netDesign', #NX must be first, netDesign is given name of costFn
                                                                 'effMetrics','resMetrics', 'os' ]), globals=tuple())
            jobs.append((netIdx,job))
        for netIdx,job in enumerate(jobs):
            ret = job[1]()
            if type(ret) != type(None):
                if self.outputState('debug'):
                    #print type(ret)
                    #print ret
                    print 'Job returned successfully.'
                    print 'ret=%s'%str(ret)
                self.netReports[netIdx]['cost'] = ret
            else:
                print type(ret)
                print 'baseDesign.evaluate(), job_number: %d'%job[0]
                print 'configParams: \n' + str(self.configParams)
                print 'costParams:   \n' + str(self.costParams)
                print 'fixedParams:  \n' + str(self.fixedParams)
                raise RuntimeError, 'Error: cost evaluation job failed!'

        if len(self.netReports) == 0:
            print 'WARNING: no networks generated!'
            return 0.0
        
        self.fitnessStats()
        return self.cost


    def evaluateHelper(self, c):
    #this creates a primitive lambda fn.  cannot be placed in a loop!
        return lambda : c

    def fitnessStats(self):
        if self.netReports == []: 
            self.evaluate()
        
        if self.fitness == None or self.cost == np.inf: 
            fitnesses    = [netReport['cost']['fitness']    for netReport in self.netReports]
            resiliences  = [netReport['cost']['resilience'] for netReport in self.netReports]
            efficiencies = [netReport['cost']['efficiency'] for netReport in self.netReports]

            self.fitness =  {'fitness':np.average(fitnesses),
                                'resilience':np.average(resiliences),
                                'efficiency':np.average(efficiencies)}
            avg = np.average(fitnesses)
            if len(self.netReports) > 1:
                cvar      = abs(np.std(fitnesses)/avg)
                self.cost = avg
                if cvar > 0.2:
                    print 'Warning: The coefficient of variation of cost = %.4G.  Perhaps increase size of sample networks?'%cvar
            else:
                self.cost = self.fitness['fitness']
            if self.outputState('debug'): 
                print self.cost
        return self.fitness

    def fitnessVector(self):
        if self.fitness == None: 
            self.fitnessStats()

        return np.array([self.fitness['resilience'], self.fitness['efficiency']])


    def mutate(self, rate=1.0):
        #override me!
        raise Exception, 'call of abstract method'


    def outputState(self, condition):
        '''
          checks to see if condition (e.g. 'debug') is in outputLevel
        '''
        if 'outputLevel' in self.fixedParams:
            return condition in self.fixedParams['outputLevel']
        else:
            return True

    def purge(self):
        self.cost = np.inf
        self.fitness = None
        self.nets = []
        self.netReports = []

    def setParam(self, paramName, newVal):
        if paramName == 'costFn':
            if not hasattr(newVal, '__call__'):
                raise ValueError, 'Parameter costFn must be a function'
            self.costParams['costFn'] = newVal
            self.purge() 
        elif paramName == 'g':
            if newVal < 0.:
                raise ValueError, 'Parameter g must >= 0.'
            self.costParams['g'] = newVal
            self.purge()
        elif paramName == 'nn':
            if newVal < 0.:
                raise ValueError, 'Number of nodes must >= 0'
            self.fixedParams['nn'] = int(round(newVal))
            self.purge()
        elif paramName == 'numNets':
            if newVal < 0.:
                raise ValueError, 'Number of nets must >= 0'
            self.fixedParams['numNets'] = int(round(newVal))
            self.purge()
        elif paramName == 'outputLevel':
            if not type(newVal) is set:
                raise ValueError, 'Parameter outputLevel must be a set'
            self.fixedParams['outputLevel'] = newVal
            self.purge() #optional
        elif paramName == 'r':
            if newVal < 0. or newVal > 1.:
                raise ValueError, 'Parameter r must lie in [0,1]'
            self.costParams['r'] = newVal
            self.purge()
        elif paramName == 'tau':
            if newVal < 0. or newVal > 1.:
                raise ValueError, 'Parameter tau must lie in [0,1]'
            self.costParams['tau'] = newVal
        else:
            raise ValueError, 'Unknown parameter: ' + paramName

    def vectorEvaluate(self, data_vector):
        #override me!
        raise Exception, 'call of abstract method'
        #return self.evaluate()

    def vectorSave(self):
        #override me!
        raise Exception, 'call of abstract method'
        #data_vector = [npr.rand(),npr.rand()]
        #return data_vector

class constantDesign(baseDesign):
    #nn = number of nodes
    #k  = mean nodes per clique

    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        if 'G' not in self.fixedParams:
            self.fixedParams['G'] = nx.Graph()
            self.fixedParams['nn'] = self.fixedParams['G'].number_of_nodes()
            print 'WARNING: initialied design to an empty graph'

    def buildNet(self):
        return self.fixedParams['G'] 

    def mutate(self):
    	return

    def setParam(self, paramName, newVal):
        if paramName == 'G':
            self.fixedParams['G']  = newVal
            self.fixedParams['nn'] = newVal.number_of_nodes()
        elif paramName == 'dummy':  #needed for grid search
            pass
        else:
            baseDesign.setParam(self, paramName, newVal)

    def vectorEvaluate(self, data_vector):
        #does nothing with the data_vector
        numTrials = 10
        for trial in xrange(numTrials):
            try:
                self.evaluate()
                break
            except Exception, inst:
                print 'Error has occurred during evaluation.'
                print 'Resubmitting.  trial %d of %d times.'%(trial+1,numTrials)
        return -self.cost

    def vectorSave(self):
        raise ValueErrror, 'vectorSave cannot be called on constDesign b/c it\'s not parametrized'


if __name__ == '__main__':
    pass
