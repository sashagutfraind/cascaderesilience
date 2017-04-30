'''
Search for solutions using
* Grid search
* Nelder-Mead

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


Warning:
    The algorithm might not too well with discrete parameters...

TODO:
    1. 'cost' should be renamed 'fitness'
'''

import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
import networkx as nx
import numpy as np
import numpy.random as npr
import random
import scipy
import scipy.optimize
import time, os, sys, pdb

def gridSearch(params):
    '''
    #wishlist: 
      1. deal with case 3D and higher meshes
      2. refinement - secondary importance
      3. parallelization
    '''
    #advantage over nelder-mead is that it's possible to parallelize at the level of the grid, rather than at the individual configuration
    #as a result, there is no need for concurrency at the thread level
    design          = params['design']
    saParams        = params.get('searchAlgParams', {})

    outputLevel     = saParams.get('outputLevel',       set())
    simulatedSearch = saParams.get('simulatedSearch',   False)
    paramData       = saParams.get('paramData', None)  #[(name1, range1), (name2, range2), ..]

    if paramData == None:
        raise ValueError, 'Empty search problem!'

    numVars = len(paramData)
    Xnames  = numVars*[None]
    Xranges = numVars*[None]
    
    if numVars == 0:
        raise ValueError, 'Empty search problem!'
    
    problemSize = 1
    for i,data in enumerate(paramData):
        Xnames[i]  = paramData[i][0]
        Xranges[i] = paramData[i][1]
        problemSize *= len(Xranges[i])
    print 'Number of grid points: %d'%problemSize 
    
        #Xmesh, Ymesh = np.meshgrid(Xs, Ys)

    Zshape = tuple([len(rng) for rng in Xranges])
    Zs = np.zeros(shape=Zshape, dtype=np.double)
    
    if numVars > 3:
        raise ValueError, 'No support for > 3D grid search'
 
    Zfull = {}
    paramVals  = {}
    badParamVals  = {}

    if 'normal' in outputLevel:
        print 'Starting Grid Search Algorithm:'
    
    best       = design(costParams=params['costParams'], fixedParams=params['fixedParams'])
    best.cost  = -np.inf
    
    curSoln   = design(costParams=params['costParams'], fixedParams=params['fixedParams'])
    extraData = {'Zs': Zs, 'Zfull':Zfull, 'paramVals':paramVals, 'badParamVals':badParamVals}

    try:
        if simulatedSearch:
            soln = curSoln.vectorSave()
            curSoln.vectorEvaluate(soln)
            best = curSoln
        else:
            if numVars == 1:
                for i, x0 in enumerate(Xranges[0]):
                    curSoln.setParam(Xnames[0], x0)
                    try:
                        ptVal    = curSoln.evaluate()
                    except Exception, e:
                        if 'debug' in outputLevel:
                            print 'construction/evaluation failed:'
                            print e
                        badParamVals[i] = {Xnames[0]: x0}
                        continue
                    Zs[i]    = ptVal  #note: Zs is indexed by i, not by x (likewise below)
                    Zfull[i] = curSoln.netReports
                    paramVals[i] = {Xnames[0]: x0}

                    if curSoln.cost > best.cost:
                        best = design(curSoln)
            elif numVars == 2:
                for i, x0 in enumerate(Xranges[0]):
                    curSoln.setParam(Xnames[0], x0)
                    for j, x1 in enumerate(Xranges[1]):
                        curSoln.setParam(Xnames[1], x1)
                        try:
                            ptVal    = curSoln.evaluate()
                        except Exception, e:
                            badParamVals[(i,j)] = {Xnames[0]: x0, Xnames[1]: x1}
                            if 'debug' in outputLevel:
                                print 'construction/evaluation failed:'
                                print e
                            continue
                        Zs[i,j]     = ptVal
                        Zfull[(i,j)]= curSoln.netReports
                        paramVals[(i,j)] = {Xnames[0]: x0, Xnames[1]: x1}

                        if curSoln.cost > best.cost:
                            best = design(curSoln)

            elif numVars == 3:  
                for i, x0 in enumerate(Xranges[0]):
                    curSoln.setParam(Xnames[0], x0)
                    for j, x1 in enumerate(Xranges[1]):
                        curSoln.setParam(Xnames[1], x1)
                        for k, x2 in enumerate(Xranges[2]):
                            curSoln.setParam(Xnames[2], x2)
                            try:
                                ptVal    = curSoln.evaluate()
                            except Exception, e:
                                if 'debug' in outputLevel:
                                    print 'construction/evaluation failed:'
                                    print e
                                badParamVals[(i,j)] = {Xnames[0]: x0, Xnames[1]: x1}
                                continue
                            Zs[i,j,k]     = ptVal
                            Zfull[(i,j,k)]= curSoln.netReports
                            paramVals[(i,j,k)] = {Xnames[0]: x0, Xnames[1]: x1, Xnames[2]: x2}

                            if curSoln.cost > best.cost:
                                best = design(curSoln)

    #conceivably, this search code could be made recursive, 
    #so for numVars >= 2, we just call outselves and integrate the results

    except RuntimeError, e:
        print e
        print 'Error: cannot complete run!'
        print 'Likely b/c network couldn\'t be completed.'
        curSoln.cost = -np.inf
        raise

    if len(Zfull) == 0:
        raise ValueError, 'Grid search failed: no configuration could be constructed/evaluated. \nXnames:\n%s\nXranges:\n%s\n'%(str(Xnames),str(Xranges))

    if 'normal' in outputLevel:
        print 'Optimal set: ', best
    sys.stdout.flush()
    
    if 'extraData' in outputLevel:
        return {'optConfig':best, 'extraData':extraData}
    else:
        return best


def nelderMead(params):
    design      = params['design']
    saParams    = params.get('searchAlgParams', {})
    runs        = saParams.get('runs',        10)
    simulatedSearch = saParams.get('simulatedSearch',   False)
    
    verbose     = saParams.get('verbose',     True)
    verbose2    = saParams.get('verbose2',    True)
    
    if verbose:
        print 'Starting Nelder-Mead Algorithm:'
        print 'Runs = %d, ' %(runs, )
    
    best      = design(costParams=params['costParams'], fixedParams=params['fixedParams'])
    best.cost = -np.inf
    
    curSoln = design(costParams=params['costParams'], fixedParams=params['fixedParams'])
    iterationData = []
    for r in xrange(runs):
        if verbose:
            print 'Run #%0.3d,'%r
        sys.stdout.flush()
        curSoln = randomizeInitialSoln(curSoln)

        try:
            if simulatedSearch:
                soln = curSoln.vectorSave()
            else:
                soln = scipy.optimize.fmin(func=curSoln.vectorEvaluate, x0=curSoln.vectorSave(), disp=1, callback=nelderMead_callback)
        except RuntimeError, e:
            print e
            print 'Error: cannot complete run!'
            print 'Likely b/c network couldn\'t be completed.'
            curSoln.cost = np.inf
            break

        curSoln.vectorEvaluate(soln)
        if curSoln.cost > best.cost:
            best = design(curSoln)

        iterationData.append((curSoln.cost, best.cost))

        if verbose:
            print 'Optimal set: ', curSoln
        sys.stdout.flush()
    
    if not verbose2:
        return best
    else:
        return {'optConfig':best, 'iterationData':iterationData}

def nelderMead_callback(curSolnVector):
    pass
    #print curSolnVector
    #maybe store into a global array...


def randomizeInitialSoln(curSoln):
    for i in xrange(40):
        curSoln.mutate()

    return curSoln

def testGridSearch():
    print 'Testing grid search algorithm...'

    import netDesign 
    #d = netDesign.constantDesign()

    class specialDesign(netDesign.baseDesign):
        def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
            netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
            
        def buildNet(self):
            nn = netDesign.baseDesign.getParam(self, 'nn')
            return nx.generators.empty_graph(nn)
            #WARNING: this test violates the convention that 'nn' is a fixedParam - cannot be configured
            #we need this since costFn only takes the graph structure, so the param must be encoded in the topology

        def evaluate(self, configParams=None, nets=None):
            return netDesign.baseDesign.evaluate(self, configParams=configParams, nets=nets)
            

    def testCostFn(G, costParamsInput, fullData=True):
        x = G.number_of_nodes() / 1000.

        fitness = x*(1-x)  #max is found at x=0.5, y=0.25

        print fitness
        r = costParamsInput.get('r', 0.5)
        if fullData:
            return {'fitness': fitness, 'resilience':fitness*r, 'efficiency':fitness*(1-r)}
        else:
            return float(fitness)

    paramSet = {}
    paramSet['design']    = specialDesign
    paramSet['costParams']= {'costFn': testCostFn}
    paramSet['fixedParams']= {'numNets':1, 'outputLevel': set(['normal'])}
    paramSet['searchAlgParams'] = {'paramData': [('nn',np.arange(0, 1000.1, 50))],
                                   'outputLevel': set(['debug']),    
                                   'simulatedSearch': False}

    optNetConfig = gridSearch(paramSet)
    assert optNetConfig.cost == 0.25
    assert optNetConfig.getParam('nn') == 500

    print 'Test passed!'


if __name__ == '__main__':
    testGridSearch()
