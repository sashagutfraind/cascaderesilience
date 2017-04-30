'''
Network Generation Programs - Multicell
    Each network consists of disconnected components: 
    * cavemen - some inter-clique connectivity
    * cliques
    * cycle
    * stars

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


'''
import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
#from netDesign import *
import netDesign

import networkx as nx
from networkx import *
import numpy as np
import numpy.random as npr
import time
import random

def buildCellularNetFractional(nn, k, cellMakerFn):
    #constructs a cellular graph, where if nn%k != 0, there is a fractionally-sized cell
    nn = int(nn)
    k  = int(round(k))
    G  = nx.Graph()
    while True:
        cell     = cellMakerFn(k)  
        G        = nx.operators.disjoint_union(G, cell)

        if G.number_of_nodes() + k > nn:
            break

    fractionCellSize = nn - G.number_of_nodes()
    if fractionCellSize > 0:
        cell     = cellMakerFn(fractionCellSize)  
        G        = nx.operators.disjoint_union(G, cell)

    return G


def buildCellularNetRandom(nn, k, cellMakerFn):
    #constructs a cellular graph, where the cell size is randomized around k (allows for approximation when nn%k!=0)
    nn = int(nn)
    k  = int(round(k))
    G  = nx.Graph()
    while True:
        cellSize = int(round(npr.normal(k, 0.3)))
        cell     = cellMakerFn(cellSize)  
        G        = nx.operators.disjoint_union(G, cell)

        if nn < G.number_of_nodes() + k/2.:
            break
    return G


#consists of weakly-connected clique graphs
class cavemenDesign(netDesign.baseDesign):
    #nn = number of nodes
    #k  = mean nodes per clique
    #p  = inter-clique connectedness probability

    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        #if x==None and costParams == None:
        #    self.costParams['costFn']  = defaultRingCost
        if x==None and fixedParams == None:
            self.fixedParams['nn']     = 100
        if x==None and configParams == None:
            self.configParams['k']     = 5.0
        if x==None and configParams == None:
            self.configParams['p']     = 0.5 #maybe better: 1./(nn*1.0/k)

            #it might be wiser to do a randomization
        if (self.fixedParams['nn'] % self.configParams['k']) != 0:
            #raise ValueError, 'Cannot initialize cavemen design when k does not evenly divide nn'
            print 'WARNING: initializing cavemen design that k does not evenly divide nn'

    def buildNet(self):
        G  = buildCellularNetFractional(nn=self.fixedParams['nn'], k=self.configParams['k'], cellMakerFn=nx.generators.complete_graph)

        p = self.configParams['p']
        ccs = connected_components(G)
        for ccNum, cc1 in enumerate(ccs):
            for cc2 in ccs[ccNum+1:]:
                if npr.rand() < p:
                    G.add_edge(cc1[0], cc2[0])  
                    #node 0 is not special since the cliques are symmetric

        return G

    '''testing:
    def evaluate(self):
        self.cost = (self.configParams['p'] -0.5) ** 2 + (self.configParams['k'] - 4.0) ** 2
        return self.cost
    '''

    def mutate(self, rate=1.0):
        if npr.rand() > rate: return

        newVal = self.configParams['k'] + self.fixedParams['nn']*(npr.rand() - 0.5) * 0.25
        self.setK(newVal)
        newVal = self.configParams['p'] + (npr.rand() - 0.5) * 0.1
        self.setP(newVal)

    	return

    def setK(self, newVal):
        newVal = max(newVal, 1.0)
        self.configParams['k'] = float(newVal)
        self.purge()

    def setP(self, newVal):
        newVal = max(newVal, 0.0)
        newVal = min(newVal, 1.0)
        self.configParams['p'] = newVal
        self.purge()

    def setParam(self, paramName, newVal):
        if paramName == 'k':
            self.setK(newVal)
        elif paramName == 'p':
            self.setP(newVal)
        else:
            netDesign.baseDesign.setParam(self, paramName, newVal)


    def vectorEvaluate(self, data_vector):
        self.setK(data_vector[0])
        self.setP(data_vector[0])

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
        data_vector = []
        data_vector.append( self.configParams['k'] )
        data_vector.append( np.tan(np.pi*self.configParams['p'] - np.pi/2.) )

        return data_vector



def eff_clique(n, k, g): 
    return 1./float(n-1) * float(k-1)

def eff_cycle(n, k, g):
    k = int(k)
    if k%2 != 0: #odd
        return  2./(n-1.) * sum([1./(j**g) for j in xrange(1,(k+1)/2)])
    else:        #even
        return  2./(n-1.) * sum([1./(j**g) for j in xrange(1,k/2)] + [0.5/(k/2.)**g])

def eff_star(n, k, g): 
    return 1.0/(n-1)*(1-1./k)*(2 + (k-2.)/(2.**g))
            

#consists of disconnected clique graphs
class cliqueDesign(netDesign.baseDesign):
    #nn = number of nodes
    #k  = mean nodes per clique

    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        #if x==None and costParams == None:
        #    self.costParams['costFn']  = defaultRingCost
        if x==None and fixedParams == None:
            self.fixedParams['nn']     = 100
        if x==None and configParams == None:
            self.configParams['k']     = 5
            #it might be wiser to do a randomization
        if (self.fixedParams['nn'] % self.configParams['k']) != 0:
            #raise ValueError, 'Cannot initialize clique design when k does not evenly divide nn'
            print 'WARNING: initializing clique design that k does not evenly divide nn'

    def buildNet(self):
        return buildCellularNetFractional(nn=self.fixedParams['nn'], k=self.configParams['k'], cellMakerFn=nx.generators.complete_graph)

    def evaluateAnalytic(self):
        self.netReports = evaluateAnalyticHelper(eff_fn=eff_clique, res_fn=res_clique, nn=self.fixedParams['nn'], k=self.configParams['k'], costParams=self.costParams)
        self.fitnessStats()
        return self.cost

    def mutate(self, rate=1.0):
        if npr.rand() > rate: return

        newVal = self.configParams['k'] + self.fixedParams['nn']*(npr.rand() - 0.5) * 0.25
        self.setK(newVal)

    	return

    def setK(self, newVal):
        newVal = max(newVal, 1.0)
        self.configParams['k'] = float(newVal)
        self.purge()

    def setParam(self, paramName, newVal):
        if paramName == 'k':
            self.setK(newVal)
        elif paramName == 'analytic':
            self.costParams['analytic'] = bool(newVal)
        else:
            netDesign.baseDesign.setParam(self, paramName, newVal)

    def vectorEvaluate(self, data_vector):
        self.setK(data_vector[0])

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
        data_vector = []
        data_vector.append( self.configParams['k'] )

        return data_vector



#consists of disconnected cycle graphs
class cycleDesign(netDesign.baseDesign):
    #nn = number of nodes
    #k  = mean nodes per cycle

    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        #if x==None and costParams == None:
        #    self.costParams['costFn']  = defaultRingCost
        if x==None and fixedParams == None:
            self.fixedParams['nn']     = 100
        if x==None and configParams == None:
            self.configParams['k']     = 5.0
        if (self.fixedParams['nn'] % self.configParams['k']) != 0:
            #raise ValueError, 'Cannot initialize cycle design when k does not evenly divide nn'
            print 'WARNING: initializing cycle design that k does not evenly divide nn'

    def buildNet(self):
        return buildCellularNetFractional(nn=self.fixedParams['nn'], k=self.configParams['k'], cellMakerFn=nx.generators.cycle_graph)

    def evaluateAnalytic(self):
        self.netReports = evaluateAnalyticHelper(eff_fn=eff_cycle, res_fn=res_cycle, nn=self.fixedParams['nn'], k=self.configParams['k'], costParams=self.costParams)
        self.fitnessStats()
        return self.cost

    def mutate(self, rate=1.0):
        if npr.rand() > rate: return

        newVal = self.configParams['k'] + self.fixedParams['nn']*(npr.rand() - 0.5) * 0.25
        self.setK(newVal)

    	return

    def setK(self, newVal):
        newVal = max(newVal, 1.0)
        self.configParams['k'] = float(newVal)
        self.purge()

    def setParam(self, paramName, newVal):
        if paramName == 'k':
            self.setK(newVal)
        elif paramName == 'analytic':
            self.costParams['analytic'] = bool(newVal)
        else:
            netDesign.baseDesign.setParam(self, paramName, newVal)


    def vectorEvaluate(self, data_vector):
        self.setK(data_vector[0])

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
        data_vector = []
        data_vector.append( self.configParams['k'] )

        return data_vector

#like cavemen, but with stars connected randomly
class connectedStarsDesign(netDesign.baseDesign):
    #nn = number of nodes
    #k  = mean nodes per clique
    #p  = inter-clique connectedness probability

    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        #if x==None and costParams == None:
        #    self.costParams['costFn']  = defaultRingCost
        if x==None and fixedParams == None:
            self.fixedParams['nn']     = 100
        if x==None and configParams == None:
            self.configParams['k']     = 5.0
        if x==None and configParams == None:
            self.configParams['p']     = 0.5

            #it might be wiser to do a randomization
        #if (self.fixedParams['nn'] % self.configParams['k']) != 0:
        #    #raise ValueError, 'Cannot initialize cavemen design when k does not evenly divide nn'
        #    print 'WARNING: initializing connected stars design that k does not evenly divide nn'

    def buildNet(self):
        star_maker_fn = lambda k: nx.generators.star_graph(k-1) #size is cellSize when including the center 
        G = buildCellularNetFractional(nn=self.fixedParams['nn'], k=self.configParams['k'], cellMakerFn=star_maker_fn)

        p = self.configParams['p']
        leaders = [node for node in G.nodes() if nx.degree(G, node) > 1]
        if leaders == []:  #occurs if k==1
            leaders = G.nodes()
        for leaderNum, leader1 in enumerate(leaders):
            for leader2 in leaders[leaderNum+1:]:
                if npr.rand() < p:
                    G.add_edge(leader1, leader2)  

        return G

    def mutate(self, rate=1.0):
        if npr.rand() > rate: return

        newVal = self.configParams['k'] + self.fixedParams['nn']*(npr.rand() - 0.5) * 0.25
        self.setK(newVal)
        newVal = self.configParams['p'] + (npr.rand() - 0.5) * 0.1
        self.setP(newVal)

    	return

    def setK(self, newVal):
        newVal = max(newVal, 1.0)
        self.configParams['k'] = float(newVal)
        self.purge()

    def setP(self, newVal):
        newVal = max(newVal, 0.0)
        newVal = min(newVal, 1.0)
        self.configParams['p'] = newVal
        self.purge()

    def setParam(self, paramName, newVal):
        if paramName == 'k':
            self.setK(newVal)
        elif paramName == 'p':
            self.setP(newVal)
        else:
            netDesign.baseDesign.setParam(self, paramName, newVal)


    def vectorEvaluate(self, data_vector):
        self.setK(data_vector[0])
        self.setP(data_vector[0])

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
        data_vector = []
        data_vector.append( self.configParams['k'] )
        data_vector.append( np.tan(np.pi*self.configParams['p'] - np.pi/2.) )

        return data_vector

def evaluateAnalyticHelper(eff_fn, res_fn, nn, k, costParams):
    if 'g' not in costParams:
        g = 1.0
        print 'Warning: using default g (%f)!'%g
    else:
        g = costParams['g']
    if 'tau' not in costParams:
        tau = 0.5
        print 'Warning: using default tau (%f)!'%tau
    else:
        tau = costParams['tau']
    if 'r' not in costParams:
        r = 0.5
        print 'Warning: using default r (%f)!'%r
    else:
        r = costParams['r']

    if nn > 1 and k > 1: 
        efficiency = eff_fn(n=nn, k=k, g=g)
        resilience = res_fn(n=nn, k=k, tau=tau)
    else:
        efficiency = 0.0
        resilience = 1.0

    fitness   = r*resilience + (1.0-r)*efficiency
    netReports = [{'cost':{'resilience':resilience, 'efficiency':efficiency, 'fitness':fitness}}]

    return netReports


def res_clique(n, k, tau):
    #based on DRAIEF, GANESH, AND MASSOULIE. J.Appl.Prob. 2008.
    import scipy.optimize
    c     = tau*(k-1)
    if c < 1:
        #wishlist: this bound on resilience is WAY too low
        return 1. - 1./(n-1.) * min(k-1, 1./(1.-c) - 1) #this is only a bound!
    elif c==1.:
        print 'WARNING: No known solution for c=1!'
        if k > 2:
            print 'Computing an upper bound...'
            return 1. - 1./(n-1.) *  (1./(1-tau*(k-2)) - 1 )
        else:
            return 1. - 1./(n-1.) *  (1.0)
            #not in the paper, but a sensible estimate
    else: #c>1
        fnc   = lambda gamma: gamma + np.exp(-c*gamma) - 1. 
        gamma = scipy.optimize.brentq(f=fnc, a=1E-10, b=1.)  #1E-10 is to avoid soln of 0.
        return  1. - 1./(n-1.) * max(1E-10,gamma) * (k-1.)

def res_cycle(n, k, tau): 
    if tau < 1.0:
        return 1-1.0/(n-1)*(2*tau*(1-(tau**(k-1.)))/(1.-tau) - (k-1.)*(tau**k) )
    else:
        return 1-(k-1.)/(n-1.)

def res_star(n, k, tau): 
    return 1-tau*(k-1.)/(n-1.)/k * (2. + tau*(k-2.))


#consists of disconnected star graphs
class starDesign(netDesign.baseDesign):
    #nn = number of nodes
    #k  = mean nodes per star

    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        #if x==None and costParams == None:
        #    self.costParams['costFn']  = defaultRingCost
        if x==None and fixedParams == None:
            self.fixedParams['nn']     = 100
        if x==None and configParams == None:
            self.configParams['k']     = 5.0
        if (self.fixedParams['nn'] % self.configParams['k']) != 0:
            #raise ValueError, 'Cannot initialize star design when k does not evenly divide nn'
            print 'WARNING: initializing star design that k does not evenly divide nn'

    def buildNet(self):
        star_maker_fn = lambda k: nx.generators.star_graph(k-1) #size is cellSize when including the center 
        return buildCellularNetFractional(nn=self.fixedParams['nn'], k=self.configParams['k'], cellMakerFn=star_maker_fn)


    def evaluateAnalytic(self):
        self.purge()
        self.netReports = evaluateAnalyticHelper(eff_fn=eff_star, res_fn=res_star, nn=self.fixedParams['nn'], k=self.configParams['k'], costParams=self.costParams)
        self.fitnessStats()
        return self.cost

    def mutate(self, rate=1.0):
        if npr.rand() > rate: return

        newVal = self.configParams['k'] + self.fixedParams['nn']*(npr.rand() - 0.5) * 0.25
        self.setK(newVal)

    	return

    def setK(self, newVal):
        newVal = max(newVal, 1.0)  #k==1 is lower bound, and when generating cells we make star(k-1)
        self.configParams['k'] = float(newVal)
        self.purge()

    def setParam(self, paramName, newVal):
        if paramName == 'k':
            self.setK(newVal)
        elif paramName == 'analytic':
            self.costParams['analytic'] = bool(newVal)
        else:
            netDesign.baseDesign.setParam(self, paramName, newVal)

    def vectorEvaluate(self, data_vector):
        self.setK(data_vector[0])

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
        data_vector = []
        data_vector.append( self.configParams['k'] )

        return data_vector

def run_fitness_k():
    #plots fitness vs. k: used for analyzing the effect of fractional cells when nn%k!=0
    nn     = 60  
    kvals  = range(1,nn+1)
    tau    = 0.8

    fitnesses    = []
    resiliences  = []
    efficiencies = []

    des = starDesign()
    des.costParams['tau'] = tau

    des.fixedParams['nn']  = nn 
    des.fixedParams['numNets']  = 1
    des.fixedParams['outputLevel'] = set()
    for k in kvals:
        des.setK(k)
        des.evaluate()

        fitnesses    += [np.average([report['cost']['fitness'] for report in des.netReports])]
        resiliences  += [np.average([report['cost']['resilience'] for report in des.netReports])]
        efficiencies += [np.average([report['cost']['efficiency'] for report in des.netReports])]

    print 'k values: '
    print kvals
    print 'fitnesses'
    print fitnesses
    print 'resiliences'
    print resiliences
    print 'efficiencies'
    print efficiencies

    import pylab
    pylab.figure()
    pylab.title('tau=%.3f'%tau)
    pylab.hold(True)
    pylab.plot(kvals, fitnesses, label='fitnesses')
    pylab.plot(kvals, resiliences, label='resiliences')
    pylab.plot(kvals, efficiencies, label='efficiencies')
    pylab.legend(loc='best')

    pylab.show()


def test_execution():
    design = cavemenDesign()
    for i in xrange(15):
        design.mutate()
    design.evaluate()
    print design

    design.buildNets()
    design.buildNets()
    

    designC1 = cliqueDesign()
    for i in xrange(15):
        designC1.mutate()
    designC1.evaluate()
    print designC1

    designC2 = cycleDesign()
    for i in xrange(15):
        designC2.mutate()
    designC2.evaluate()
    print designC2

    designS = starDesign()
    for i in xrange(15):
        designS.mutate()
    designS.evaluate()
    print designS

if __name__ == '__main__':
    test_execution()

