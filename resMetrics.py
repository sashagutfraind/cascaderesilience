'''
Computes the resilience of a graph to an epidemic

    Copyright (C) 2007-2013 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


Notes:
-----

includes:
    * various implementations of various measures;  invokes resMetricsX.py
    * adaptive sampling is used to determine the number of runs to do, in order to ensure sufficient number of runs
    * tests

no guarantee of support for directed graphs:
    the definition of neighbors might need to change


'''

import sys, os
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
#sys.path.append('./extensions')

import timeit, pdb
import numpy as np
import numpy.random as npr
import networkx as nx

def extentSIRSimLAdaptive(G, params={}):
#list-dict implementation of the contagion. much faster than a queue or a list-list
#supports 2FreeSIR using transmission model parameter
    revisedParams = params.copy()
    nn              = G.number_of_nodes()
    tau		        = revisedParams.get('tau', 0.5)
    minReplications = revisedParams.get('minReplications', max(nn/10, 40))
    startNodes      = revisedParams.get('startNodes', range(nn))
    tolerance       = revisedParams.get('tolerance', 0.05)
    transmissionModel = revisedParams.get('transmissionModel', 'SIR')  #KFleeSIR
    vaccinated      = revisedParams.get('vaccinated', [])
    weighted        = revisedParams.get('weighted', False)
    if 'minReplications' not in revisedParams:
        print 'Warning: using detault minReplications=%d!'%minReplications
    if 'tau' not in revisedParams:
        revisedParams['tau'] = 0.5
        print 'Warning: using detault tau=%.3f!'%revisedParams['tau']
    if 'tolerance' not in revisedParams:
        print 'Warning: using detault tolerance=%.3f!'%tolerance
    if transmissionModel == 'SIR':
        extentFunc = extentSIRSim_helper
    else:
        raise ValueError, 'Unsupported epidemic: '+transmissionModel

    if G.nodes() == []: return 0., 0.

    acceptable = lambda outcomes, norm, tolerance: 1.96*np.std(outcomes) / norm / np.sqrt(len(outcomes)) < tolerance

    extents = []
    trial   = 1
    while (trial <= minReplications or not acceptable(extents, float(nn), tolerance)): # and trial < maxReplications:
        #too slow: startNode = G.nodes()[npr.randint(nn)]
        startNode = startNodes[npr.randint(len(startNodes))]
        extents  += [extentFunc(G=G, startN=startNode, params=revisedParams, vaccinated=vaccinated, weighted=weighted)]
        trial    += 1

    #print 'Num trials: %d'%trial 
    return np.average(extents), np.std(extents)

def extentSIRSim_helper(G, startN, params, vaccinated, weighted=False, verbose=False):
    tau = params['tau']
    uninfectable = dict.fromkeys(vaccinated, 'vaccinated')
    uninfectable[startN] = 'infected'
   
    infected    = G.number_of_nodes()*[None]
    infected[0] = startN
    head        = 0
    tail        = 1
    if weighted:
        while head < tail:
            curNode = infected[head]
            head   += 1
            for nb in G.neighbors(curNode):
                wt = G.get_edge_data(curNode, nb)['weight']
                if not (nb in uninfectable) and tau/wt > npr.rand():
                    uninfectable[nb] = 'infected'
                    infected[tail]   = nb
                    tail            += 1
    else:
        while head < tail:
            curNode = infected[head]
            head   += 1
            for nb in G.neighbors(curNode):
                if not (nb in uninfectable) and tau > npr.rand():
                    uninfectable[nb] = 1
                    infected[tail]   = nb
                    tail            += 1
    if verbose:
        return len(uninfectable), uninfectable
    else:
        return len(uninfectable)

    #attempts at optimizations do not seem to produce results:
    #replacing infected/uninfectable with numpy arrays, replacing uninfectable with set


def extentSIRSimX(G, params):
    import resMetricsX
    return resMetricsX.extentSIR_adaptive(G=G, params=params)

def profileCythonSim():
    params = {'tau':0.1, 'tolerance':0.05, 'minReplications':100}

    print 'Python implementation:'
    G = nx.generators.erdos_renyi_graph(100,0.3)
    t = timeit.default_timer() 
    print extentSIRSimLAdaptive(G, params=params)
    print 'Time:' + str(timeit.default_timer()-t)

    print 'Cython implementation:'
    t = timeit.default_timer() 
    try:
      print extentSIRSimX(G, params=params)
    except Exception, e:
      print e
    print 'Time:' + str(timeit.default_timer()-t)


def testCython():
    import pylab

    meansCyt = []
    sdevsCyt = []
    meansNum = []
    sdevsNum = []
    params = {'minReplications':5, 'tolerance':0.1, 'tau':0.1}
    for graphNum in xrange(100):
        G = nx.generators.erdos_renyi_graph(100,0.2)
        for trial in xrange(10):
            m, s  = extentSIRSimX(G, params)
            meansCyt.append(m)
            sdevsCyt.append(s)
            m, s  = extentSIRSimLAdaptive(G, params)
            meansNum.append(m)
            sdevsNum.append(s)
    pylab.show()
    figMC = pylab.figure()
    pylab.hold(True)
    bins = range(0,101,10)
    countsMeanCyt = pylab.hist(meansCyt, bins=bins, histtype='step', label='Cyt')[0]
    countsMeanNum = pylab.hist(meansNum, bins=bins, histtype='step', label='Num')[0]
    pylab.title('means') 

    print
    print 'Means'
    meansDiff = sum(np.abs(countsMeanCyt-countsMeanNum))
    relError = 0.  #aside: rel error increases when tau decreases -> smaller extents
    for i, mn in enumerate(meansCyt):
        if mn > 0.:
            relError += (mn-meansNum[i])/mn
    print 'Absolute difference (always>0): %f'%meansDiff
    print 'Avg relative error: %f'%(relError/len(meansCyt))
    #print 'Net      difference: %f'%sum(      (countsMeanCyt-countsMeanNum))
    if (meansDiff / float(len(bins)*(graphNum+1))) < 1 and relError < 0.2:
        print 'Test passed!'
        print 'Visually inspect: the histograms should broadly overlap'
    else:
        print 'Test failed: more than one misplaced count per graph per bin!'


    print
    print 'SDevs'
    pylab.show()
    figMC = pylab.figure()
    pylab.hold(True)
    bins = range(0,101,10)
    countsStdCyt = pylab.hist(sdevsCyt, bins=bins, histtype='step', label='Cyt')[0]
    countsStdNum = pylab.hist(sdevsNum, bins=bins, histtype='step', label='Num')[0]
    pylab.title('sdevs') 

    stdevsDiff = sum(np.abs(countsStdCyt-countsStdNum))
    relError = 0.
    for i, sd in enumerate(sdevsCyt):
        if sd > 0.:
            relError += (sd-sdevsNum[i])/mn
    print 'Absolute difference: %f'%stdevsDiff
    print 'Avg of relative error: %f'%(relError/len(sdevsCyt))



if __name__ == '__main__':
    print 'Running Tests...\n'
    profileCythonSim()
    testCython()
    #testAnalytic()
