'''
Simulations of SI on cliques

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 
 
#TODO
0. Plug SI soln back into the ODE
1. Continuous time SI
    c) compute mean time to becoming infected, compare to analytic
2. Continuous time SIR
    a) compute average extent
    b) plot avg extent vs. p, mu

'''

import os, sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
import networkx as nx
#from networkx.centrality import *
import numpy as np
import numpy.random as npr
import time
import random
import pylab

import resMetrics
import effMetrics
import graphStats

analyticSI = lambda nns, beta, t: (1.-np.exp(-beta*t*(nns-1)))/(1.+(nns-2)*np.exp(-beta*t*(nns-1)))

def simulateSIR(params, delT=0.01, normalizedExtents=False):
    nn    = params.get('nn', 5)
    beta  = params.get('beta',  0.4) * delT
    g     = params.get('g', 0.1) * delT
    repetitions = params.get('repetitions', 50)

    extents = [] 
    iSeries = []
    rSeries = []

    for run in xrange(repetitions):
        startNode = 0
            
        infected  = {startNode:1} #I
        recovered = {} #R

        iSeries.append([0])
        rSeries.append([0])

        while infected != {}:
            newCases = {}
            newRecovered = {}
            for curNode in infected:
                for nb in range(nn):
                    if not (nb in infected or nb in recovered) and beta > npr.rand():
                        newCases[nb] = 1
                        #print str(curNode) + ' infected '+str(nb)
                if g > npr.rand():
                    newRecovered[curNode] = 1
                    #print 'Cured '+str(curNode)
            infected .update(newCases)
            recovered.update(newRecovered)
            for node in newRecovered:
                infected.pop(node)

            iSeries[run].append((len(infected)-1.)/(nn-1.))
            rSeries[run].append((len(recovered)-1.)/(nn-1.))
        
        if normalizedExtents:
            extents.append((len(recovered)-1.)/(nn-1.))
        else:
            extents.append(len(recovered))
            
    return extents, iSeries, rSeries


def simulateSI(params, delT=0.01, normalizedExtents=False):
    nn   = params.get('nn', 5)
    beta = params.get('beta', 0.4) * delT
    repetitions = params.get('repetitions', 50)

    extents = [] 
    iSeries = []

    timesAt50  = np.zeros(repetitions, dtype=np.double) #50% at I
    timesAt90  = np.zeros(repetitions, dtype=np.double) #90% at I

    for run in xrange(repetitions):
        zeit = 0.
        startNode = 0
            
        infected  = {startNode:1} #I
    
        iSeries.append([0.])
        
        while infected != {}:
            newCases = {}
            for curNode in infected:
                for nb in range(nn):
                    if not (nb in infected) and beta > npr.rand():
                        newCases[nb] = 1
                        #print str(curNode) + ' infected '+str(nb)
            infected.update(newCases)

            p = (len(infected)-1.)/(nn-1.)
            iSeries[run].append(p)
            if timesAt50[run] == 0. and p > 0.5:
                timesAt50[run] = zeit
            if timesAt90[run] == 0. and p > 0.9:
                timesAt90[run] = zeit
            if abs(p-1.) < .1/nn:
                break

            zeit += delT
        
        if normalizedExtents:
            extents.append((len(infected)-1.)/(nn-1.))
        else:
            extents.append(len(infected))
            
    return extents, iSeries, timesAt50, timesAt90


def testAnalyticSI(nn=100, beta=0.1, repetitions=50):
    delT   = 0.01
    params = {'nn':nn, 'beta':beta, 'repetitions':repetitions}

    extents, iSeries, timesAt50, timesAt90 = simulateSI(params, delT,normalizedExtents=True)

    pylab.figure(0)
    pylab.hold(True)
    simLengths = []
    for run in range(len(iSeries)):
        simLengths.append(len(iSeries[run]))

    tmax  = max(simLengths)
    times = delT*np.arange(0, tmax, delT)
    analyticSIseries = analyticSI(nn,beta,times)
    pylab.plot(times, analyticSIseries, 'r-', linewidth=4.)

    for run in range(len(iSeries)):
        runTimes = delT*np.arange(0, simLengths[run])
        pylab.plot(runTimes, iSeries[run])

    pylab.show()


def plotSIRsims(nn=100, beta=0.1, g=0.7, repetitions=50):
    delT   = 0.01
    params = {'nn':nn, 'beta':beta, 'g':g, 'repetitions':repetitions}

    extents, iSeries, rSeries = simulateSIR(params, delT, normalizedExtents=True)

    pylab.figure(0)
    pylab.hold(True)
    simLengths = []
    for run in range(len(iSeries)):
        simLengths.append(len(iSeries[run]))

    #tmax  = max(simLengths)
    #times = delT*np.arange(0, tmax, delT)
    #analyticSIseries = analyticSI(nn,beta,times)
    #pylab.plot(times, analyticSIseries, 'r-', linewidth=4.)

    for run in range(len(iSeries)):
        runTimes = delT*np.arange(0, simLengths[run])
        pylab.plot(runTimes, iSeries[run])

    pylab.show()


def collectSIRsims(cellSizes=range(2,21), beta=0.1, g=1.0, repetitions=10, doPlot=True):
    delT   = 0.01
    params = {'beta':beta, 'g':g, 'repetitions':repetitions}

    extents = []
    for i,nn in enumerate(cellSizes):
        params['nn'] = nn
        extentsNN, iSeries, rSeries = simulateSIR(params, delT, normalizedExtents=doPlot)
        extents.append(extentsNN)

    if doPlot:
        pylab.boxplot(np.array(extents).transpose())
        pylab.xticks(range(1,len(cellSizes)), cellSizes)
        pylab.xlabel('Cell Size')
        pylab.ylabel('Fraction captured (beta=%.2f)'%beta)
        #pylab.ylabel('Captured (beta=%.2f)'%beta)
        pylab.savefig('output' + os.sep + 'cellSize' + str(time.localtime()[3:-1]) + '.eps')
        
    return extents

def plotSIRtradeoff(repetitions=10, betas=[0.1, 0.2, 0.3, 0.5, 0.9]):
    
    cellSizes=np.arange(2,21,dtype=np.double)
    pylab.figure()
    pylab.hold(True)
    pylab.xlabel('Cell Size')
    pylab.ylabel('Fraction cell surviving')
    for beta in betas:
        extents = collectSIRsims(cellSizes=cellSizes, beta=beta, repetitions=100, doPlot=False)
        means   = np.array(extents).mean(1)

        pylab.plot(cellSizes, (cellSizes-means)/cellSizes, label='beta=%.2f'%beta)
    pylab.legend(loc='center right')

'''
def computeNets(pvals = [], pathDir = 'output'):
    repeats = 5
    design  = cliqueDesign

    for p in pvals:
        print 'The p=' + str(p) + ' case:'
        params['cost']       = netDesign.defaultCost
        params['costParams'] = {'T':p}
        
        f = open(pathDir + os.sep + 'p=%.2f_'%p + 'summary_' + 'time=' + str(time.localtime()[3:-1]) + '.txt', 'w')
        dummyNetDesign = design()
        f.write(dummyNetDesign.designParamNamesStr())
        f.write(dummyNetDesign.statsNamesStr() + '\n')#os.linesep)

        for r in range(repeats):
            optNetDesign = sa.runSA(params)
            designVals = optNetDesign.designParamValsStr()
            for i in range(len(optNetDesign.nets)):
                netData = optNetDesign.statsValsStr(i)
                f.write(designVals + netData + '\n')#os.linesep)
            
        f.close()
        print 'The optimal design is:'
        print optNetDesign
        print 
'''

if __name__ == '__main__':
     testAnalyticSI(repetitions=1)

