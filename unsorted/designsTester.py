'''
Testers for all designs (general functionality, including metrics)

    Copyright (C) 2007-2009 by  
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
#for running scripts silently:
if 'matplotlib.backends' not in sys.modules: matplotlib.use('pdf')
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('PS')
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import networkx as nx
#import pylab
import time, os, sys
import cPickle

import netDesign
#import netDesignUnrooted
#import netDesignCentrality
import netDesignGnp
import netDesignMulticell

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
#timeNow = lambda : str(time.localtime()[:-6]) + '_' + str(time.localtime()[-6:-3])

def creationTest(designClass, params, expCost=None):
    print('Building design %s...'%str(designClass))
    d = designClass()
    for paramName in params:
        d.setParam(paramName, params[paramName])

    print('Building networks...')
    d.buildNets()
    for paramName in params:
        assert d.getParam(paramName) == params[paramName], 'Design %s failed to correctly set parameter %s'%(str(design), paramName)

    print('Evaluating networks...')
    d.evaluate()

    if expCost != None:
        assert np.abs(d.cost-expCost) < 5E-2 #5% tolerance

    return d


def getterSetterTest(designClass, params=None):
    print('Testing design %s...'%str(designClass))
    d = designClass()
    for paramName in params:
        d.setParam(paramName, params[paramName])
        assert d.cost == np.inf
        assert d.nets == []
        assert d.netReports == []

    return d



def testCliques():
    print 'Testing: Cliques design'
    designClass = netDesignMulticell.cliqueDesign

    paramSets = []
    costs     = []
    
    paramSets += [{'nn':100, 'k':10, 'r':0.7, 'tau':1.0, 'g':1, 'outputLevel':set()}]
    costs     += [0.7*(1.-10/99.) + 0.3*(9./99.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 1
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(0.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 100
    curParamSet['r'] = 1.0    
    paramSets += [curParamSet]
    costs     += [1.0*(0.) + 0.0*(0.5)]

    curParamSet = paramSets[0].copy()
    curParamSet['tau'] = 0.0
    curParamSet['k'] = 100
    curParamSet['analytic'] = False
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]

    curParamSet = paramSets[-1].copy()
    curParamSet['analytic'] = True
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]

    for i, curParamSet in enumerate(paramSets):
        getterSetterTest(designClass=designClass, params=curParamSet)
        creationTest(designClass=designClass, params=curParamSet, expCost=costs[i])

def testCycles():
    print 'Testing: Cycles design'
    designClass = netDesignMulticell.cycleDesign

    paramSets = []
    costs     = []
    
    paramSets += [{'nn':100, 'k':1, 'r':0.7, 'tau':1.0, 'g':1, 'outputLevel':set()}]
    costs     += [0.7*1. + 0.3*0.]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 100
    curParamSet['g'] = 0.
    curParamSet['analytic'] = False
    paramSets += [curParamSet]
    costs     += [0.7*0. + 0.3*1]

    curParamSet = paramSets[-1].copy()
    curParamSet['analytic'] = True
    paramSets += [curParamSet]
    costs     += [0.7*0. + 0.3*1]

    for i, curParamSet in enumerate(paramSets):
        getterSetterTest(designClass=designClass, params=curParamSet)
        creationTest(designClass=designClass, params=curParamSet, expCost=costs[i])


def testConnStars():
    print 'Testing: Connected Stars design'
    designClass = netDesignMulticell.connectedStarsDesign

    paramSets = []
    costs     = []
    
    paramSets += [{'nn':100, 'k':10, 'p':1.0, 'r':0.7, 'tau':1.0, 'g':0., 'outputLevel':set()}]
    costs     += [0.7*(0.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 1
    paramSets += [curParamSet]
    costs     += [0.7*(0.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['p'] = 0
    paramSets += [curParamSet]
    costs     += [0.7*(1.-10/99.) + 0.3*(9./99.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 100
    curParamSet['p'] = 0.0    
    paramSets += [curParamSet]
    costs     += [0.7*(0.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['tau'] = 0.0
    curParamSet['k'] = 100
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['tau'] = 0.0
    curParamSet['g'] = 300.
    curParamSet['k'] = 100
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(0.)]  #efficiency is nearly 0.

    for i, curParamSet in enumerate(paramSets):
        getterSetterTest(designClass=designClass, params=curParamSet)
        creationTest(designClass=designClass, params=curParamSet, expCost=costs[i])

def testCavemen():
    print 'Testing: Connected Cliques = Cavemen design'
    designClass = netDesignMulticell.cavemenDesign

    paramSets = []
    costs     = []
    
    paramSets += [{'nn':100, 'k':10, 'p':1.0, 'r':0.7, 'tau':1.0, 'g':0., 'outputLevel':set()}]
    costs     += [0.7*(0.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 1
    paramSets += [curParamSet]
    costs     += [0.7*(0.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['p'] = 0
    paramSets += [curParamSet]
    costs     += [0.7*(1.-10/99.) + 0.3*(9./99.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 100
    curParamSet['p'] = 0.0    
    paramSets += [curParamSet]
    costs     += [0.7*(0.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['tau'] = 0.0
    curParamSet['k'] = 100
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['tau'] = 0.0
    curParamSet['k'] = 100
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]  #efficiency is nearly 0.

    for i, curParamSet in enumerate(paramSets):
        getterSetterTest(designClass=designClass, params=curParamSet)
        creationTest(designClass=designClass, params=curParamSet, expCost=costs[i])

def testConstantDesign():
    print 'Testing: Constant design'
    designClass = netDesign.constantDesign

    paramSets = []
    costs     = []
    
    paramSets += [{'nn':100, 'G':nx.generators.complete_graph(100), 'r':0.7, 'tau':1.0, 'g':1, 'outputLevel':set()}]
    costs     += [0.7*(0.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['tau'] = 0.
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]

    for i, curParamSet in enumerate(paramSets):
        getterSetterTest(designClass=designClass, params=curParamSet)
        creationTest(designClass=designClass, params=curParamSet, expCost=costs[i])
    
def testFitnessFunction():
    def compareValues(fitnessData, fitness=None, efficiency=None, resilience=None):
        if fitness    != None: assert abs(fitnessData['fitness']-fitness) < 0.01
        if efficiency != None: assert abs(fitnessData['efficiency']-efficiency) < 0.01
        if resilience != None: assert abs(fitnessData['resilience']-resilience) < 0.01

    G           = nx.generators.star_graph(10)  #10 spokes nodes, 1 leader node
    fitnessData = netDesign.defaultSIRcost(G=G, costParamsInput={'g':100., 'tau':1.0, 'minReplications':40, 'tolerance':0.050, 'r':0.5}, fullData=True)
    compareValues(fitnessData, fitness=0.5*2*10./(11.*10.), efficiency=2*10./(11.*10.), resilience=0.0)

    G           = nx.generators.star_graph(10)
    fitnessData = netDesign.defaultSIRcost(G=G, costParamsInput={'g':0., 'tau':0.0, 'minReplications':40, 'tolerance':0.050, 'r':0.5}, fullData=True)
    compareValues(fitnessData, fitness=1.0, efficiency=1.0, resilience=1.0)

    G           = nx.operators.disjoint_union(nx.generators.star_graph(10), nx.generators.star_graph(10))
    fitnessData = netDesign.defaultSIRcost(G=G, costParamsInput={'g':100., 'tau':1.0, 'minReplications':400, 'tolerance':0.050, 'r':0.2}, fullData=True)
    compareValues(fitnessData, fitness=0.2*(1-10./21) + (1-0.2)*2*2*10./(22.*21.), efficiency=2*2*10./(22.*21.), resilience=1-10./21)

    G           = nx.generators.cycle_graph(100)
    fitnessData = netDesign.defaultSIRcost(G=G, costParamsInput={'g':0., 'tau':0.2, 'minReplications':400, 'tolerance':0.050, 'r':0.5}, fullData=True)
    compareValues(fitnessData, fitness=None, efficiency=1.0, resilience=1. - 2*0.2/(1-0.2)/99.)  #fluctuations sometimes set the resilience value out of tolerance


def testGnp():
    print 'Testing: G(n,p) design'
    designClass = netDesignGnp.gnpDesign

    paramSets = []
    costs     = []
    
    paramSets += [{'nn':100, 'p':0.0, 'r':0.7, 'tau':1.0, 'g':0., 'outputLevel':set()}]
    costs     += [0.7*(1.) + 0.3*(0.)]

    curParamSet = paramSets[0].copy()
    curParamSet['p'] = 0.2 #well after the transition
    curParamSet['g'] = 0.0
    paramSets += [curParamSet]
    costs     += [0.7*(0.) + 0.3*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['p'] = 0.3 #well after the transition
    curParamSet['g'] = 0.0
    curParamSet['tau'] = 0.0
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]

    for i, curParamSet in enumerate(paramSets):
        getterSetterTest(designClass=designClass, params=curParamSet)
        creationTest(designClass=designClass, params=curParamSet, expCost=costs[i])

def testMetricsAnalytic():
    '''Compares the metrics of efficiency and resilience as computed by the analytic formulas to the results of simulations on actual graphs'''
    #designs = [netDesignMulticell.cliqueDesign,netDesignMulticell.cycleDesign,netDesignMulticell.starDesign]
    #WARNING: clique design fails to compute resilience correctly when k==nn
    designs = [netDesignMulticell.cycleDesign,netDesignMulticell.starDesign]

    taus = [0.01, 0.1]
    for tau in taus:
        print 'tau: %f'%tau
        for design in designs:
            print 'Testing: ' + str(design)
            testMetricsAnalyticHelper(design=design, tau=tau)


def testMetricsAnalyticHelper(design, tau):
    tau = tau
    g   = 1.0
    nn  = 100

    kvals = [1, 2, 5, 10, 20, 25, 50, 100]  #divisors of nn.  otherwise analytic formulas may be inaccurate
    eff_anal = []
    res_anal = []
    eff_siml = []
    res_siml = []

    des = design()
    des.setParam('tau', tau)
    des.setParam('g', g)

    des.setParam('nn', nn)
    des.setParam('numNets',  20)
    des.setParam('outputLevel', set())
    for k in kvals:
        print 'k: %.2f'%float(k)
        des.setK(k)

        des.costParams['analytic'] = True
        des.evaluate()
        efficiencies = [report['cost']['efficiency'] for report in des.netReports]
        resiliences  = [report['cost']['resilience'] for report in des.netReports]
        eff_anal += [np.average(efficiencies)]
        res_anal += [np.average(resiliences)]

        des.costParams['analytic'] = False
        des.evaluate()
        efficiencies = [report['cost']['efficiency'] for report in des.netReports]
        resiliences  = [report['cost']['resilience'] for report in des.netReports]
        eff_siml += [np.average(efficiencies)]
        res_siml += [np.average(resiliences)]

    print 'k values: '
    print kvals
    print 'Analytic solutions'
    print 'Eff'+ str(eff_anal)
    print 'Res'+ str(res_anal)
    print 'Simulated solutions'
    print 'Eff'+ str(eff_siml)
    print 'Res'+ str(res_siml)

    for counter, k in enumerate(kvals):
        assert abs(eff_anal[counter]-eff_siml[counter]) < 0.05
        assert abs(res_anal[counter]-res_siml[counter]) < 0.05

    import pylab
    pylab.figure()
    pylab.hold(True)
    pylab.plot(kvals, eff_anal, label='analytic')
    pylab.plot(kvals, eff_siml, label='simulated')
    pylab.legend(loc='best')
    pylab.title('efficiencies')

    pylab.figure()
    pylab.hold(True)
    pylab.plot(kvals, res_anal, label='analytic')
    pylab.plot(kvals, res_siml, label='simulated')
    pylab.legend(loc='best')
    pylab.title('resiliences')
    pylab.show()


def testStars():
    print 'Testing: Stars design'
    designClass = netDesignMulticell.starDesign

    paramSets = []
    costs     = []
    
    paramSets += [{'nn':100, 'k':10, 'r':0.7, 'tau':1.0, 'g':0, 'outputLevel':set()}]
    costs     += [0.7*(1.-10/99.) + 0.3*(9./99.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 1
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(0.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 100
    curParamSet['r'] = 1.0    
    paramSets += [curParamSet]
    costs     += [1.0*(0.) + 0.0*(1.)]

    curParamSet = paramSets[0].copy()
    curParamSet['k'] = 100
    curParamSet['r'] = 0.0    
    curParamSet['g'] = 2.0    
    paramSets += [curParamSet]
    costs     += [0.0*(0.) + 1.0*(0.25)]

    curParamSet = paramSets[0].copy()
    curParamSet['tau'] = 0.0
    curParamSet['k'] = 100
    curParamSet['analytic'] = False
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]

    curParamSet = paramSets[-1].copy()
    curParamSet['analytic'] = True
    paramSets += [curParamSet]
    costs     += [0.7*(1.) + 0.3*(1.)]

    for i, curParamSet in enumerate(paramSets):
        getterSetterTest(designClass=designClass, params=curParamSet)
        creationTest(designClass=designClass, params=curParamSet, expCost=costs[i])


if __name__ == '__main__':
    testCliques()
    testCycles()
    testConnStars()
    testCavemen() 
    testConstantDesign()
    testGnp()
    testStars()

    testFitnessFunction()
    testMetricsAnalytic()

