'''
generation of growing networks

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
#for running scripts silently:
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('pdf')
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('PS')
#import matplotlib.pylab as pylab
#import matplotlib.pyplot as plt
import networkx as nx
import pylab  #FIXME
import time, os, sys
import cPickle

import resMetrics, effMetrics

def simpleSIRcost(G, costParamsInput, fullData=True):
#workaround resMetricsX compilation problem

    curWorkDir = costParamsInput.get('curWorkDir', None)
    if curWorkDir != None:
        os.chdir(curWorkDir)
    costParams = costParamsInput.copy()
    if 'tau' not in costParamsInput:
        costParams['tau'] = 0.5
        print 'Warning: using detault tau (%f)!'%costParams['tau']
    avgEpidemic, stdevEpidemic = resMetrics.extentSIRSim(G, params=costParams)
    #we could consider the stdevEpidemic: higher variance -> worst.  but how precisely?

    nn = G.number_of_nodes()
    efficiency = effMetrics.ringReach(G, costParams)*1.0/nn/(nn-1)
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


def secretInfoCost(G, costParamsInput, fullData=True):
#cost according to Lindelauf et al.
    curWorkDir = costParamsInput.get('curWorkDir', None)
    if curWorkDir != None:
        os.chdir(curWorkDir)
        
    costParams = costParamsInput.copy()

    nn = G.number_of_nodes()
    ne = G.number_of_edges()

    secrecy = float(nn**2 - nn - 2*ne) / nn**2

    distances = nx.all_pairs_shortest_path_length(G)
    totalDist = 0.
    for node in distances:
        if len(distances[node]) < nn:
            totalDist += np.inf
            break
        totalDist += reduce(lambda x,y:x+y, distances[node].values()) 

    information = float(nn*(nn-1)) / totalDist
    fitness = secrecy * information

    #print 'defaultSIRcost completed computation!'
    if fullData:
        return {'resilience':float(secrecy), 'efficiency':float(information), 'fitness':float(fitness),
                'secrecy':float(secrecy), 'information':float(information), 'mu':float(fitness)}
    else:
        return float(fitness)


def secretInfoCostTest():
    trueMu = []
    compMu = []
    numNodes = range(2,20)
    for k in numNodes:
        starK = nx.generators.star_graph(k-1)
        nn = starK.number_of_nodes()
        ne = starK.number_of_edges()

        trueMu += [float(nn**2 - nn - 2*ne)/(nn**2) * float(nn**2 - nn)/float(2*((nn-1)**2))]
        compMu += [secretInfoCost(starK,{})['mu']]

    pylab.figure()
    pylab.hold(True)
    pylab.plot(numNodes, trueMu, 'b-')
    pylab.plot(numNodes, compMu, 'r-')
    pylab.show()

def growNetwork(params):
#we could keep track of seeds in case we want to recreate the network
    fitnessFn = params.get('fitnessFn',  simpleSIRcost)
    seed      = params.get('seed', None)
    if seed != None:
        npr.seed(seed)

    G = nx.Graph()
    G.add_node(0)
    for i in xrange(1,100):
        j = npr.randint(i)
        G.add_edge(i,j)
    fitness = fitnessFn(G, {'tau':0.2, 'g':1.0, 'r':0.49, 'minReplications':40, 'tolerance':0.05})
    
    return {'G':G, 'fitness':fitness, 'seed':seed}

def plotDist():
    perfAr = []
    for netNum in xrange(1000):
        res = growNetwork({'seed':netNum})
        perfAr.append(res['fitness']['fitness'])
    
    perfAr = np.array(perfAr)
    #perfAr.sort()

    pylab.hist(perfAr, bins=20)

    pylab.show()

    seed_min = perfAr.argmin() #exploing: seed=indx
    #seed_med = np.nonzero(abs(perfAr-np.round(np.median(perfAr)))<0.1)
    seed_med = np.argsort(perfAr)[len(perfAr)/2]
    seed_max = perfAr.argmax()

    pylab.figure()
    print 'min F=%.4f'%perfAr[seed_min]
    G_min = growNetwork({'seed':seed_min})['G']
    nx.draw(G_min)

    pylab.figure()
    print 'med F=%.4f'%perfAr[seed_med]
    G_med = growNetwork({'seed':seed_med})['G']
    nx.draw(G_med)

    pylab.figure()
    print 'max F=%.4f'%perfAr[seed_max]
    G_max = growNetwork({'seed':seed_max})['G']
    nx.draw(G_max)


def plotEvolution(params={}):
    fitnessFn = params.get('fitnessFn',  simpleSIRcost)
    seed      = params.get('seed', None)
    if seed != None:
        npr.seed(seed)

    fitnesses = []
    G = nx.Graph()
    G.add_node(0)
    for i in xrange(1,200):
        j = npr.randint(i)
        G.add_edge(i,j)
        fitnesses += [fitnessFn(G, {'tau':0.0, 'g':1.0, 'r':0.49, 'minReplications':40, 'tolerance':0.05})['fitness']]  #FIXME: no cascade
    
    pylab.plot(fitnesses)
    return {'G':G, 'seed':seed}

    
def _test():
    #plotDist()
    plotEvolution()

if __name__ == '__main__':
    _test()

