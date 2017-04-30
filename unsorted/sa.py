'''
Simulated Annealing Solver

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 
 

'''
import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
from   networkx 	import *
import numpy
import numpy.random as npr
import netDesign

def getCost(X):
    return X.evaluate()

def getTemperatures(X, design, Pinit, Pfinal, params):
    samples = []
    failures = 0
    numSamples  = 10
    costs = numpy.zeros(numSamples, dtype=numpy.double)
    for i in xrange(numSamples):
        Xprime = design(X)
        for j in xrange(10):
            Xprime.mutate()
        try:
            Xprime.evaluate()
        except RuntimeError:
            print 'Failed to construct sample net.  using a cheat.'
            Xprime.cost = costs[i-1]
        costs[i] = Xprime.cost

    #we take the minimum of the samples, guaranteeing that all moves will be downwards - rejectable
    delcosts = costs - numpy.repeat(costs.min(), numSamples)
    costSigma = numpy.sqrt(numpy.var(delcosts))

    k = 3
    Tinit  = -k*costSigma / numpy.log(Pinit)
    Tfinal = -k*costSigma / numpy.log(Pfinal)

    #for small flat graphs, White's algorithm above may give 0,0, so we do:
    if delcosts.max() < 1E-10:
        print 'Warning! Temperature for simulated annealing may have been computed incorrectly!'
        alpha  = params.get('alpha',       1.5)
        Tinit  += alpha**3
        Tfinal += alpha**(-7) 

    return Tinit, Tfinal

def randomizeInitialX(X):
    for i in xrange(40):
        X.mutate()

    return X

def runIteration(X, design, T, params):
    #we need the cost to compare it with Xprime
    if X.cost == None:
       X.evaluate()

    Xprime = design(X)
    Xprime.mutate()
    Xprime.cost = getCost(Xprime)

    deltaCost = Xprime.cost - X.cost
    if deltaCost < 0 or numpy.exp(-deltaCost/T) > npr.uniform():
       return Xprime
    else:
       return X

def runSA(params):
    design      = params['design']
    saParams    = params.get('searchAlgParams', {})
    runs        = saParams.get('runs',        10)
    simulated   = saParams.get('simulated',   False)

    alpha       = saParams.get('alpha',       1.2)
    beta        = saParams.get('beta',        1.1)
    Minit       = saParams.get('Minit',       1.0)
    Pinit       = saParams.get('Pinit',       0.95)
    Pfinal      = saParams.get('Pfinal',      0.05)
    
    verbose     = saParams.get('verbose',     True)
    verbose2    = saParams.get('verbose2',    False)
    
    if verbose:
        print 'Starting Simulated Annealing Algorithm:'
        print 'Runs = %d, alpha = %.2f, beta=%.2f, Minit=%d, Pinit = %.2f, Pfinal = %.2f' %(runs, alpha, beta, Minit, Pinit, Pfinal)
    
    if alpha - 1.0 < 1E-5:
        raise ValueError, 'alpha is too low.  Must be > 1.0'

    best = design(costParams=params['costParams'], fixedParams=params['fixedParams'])
    best.cost       = np.inf
    bestDiscoveries = 0
    
    #we compute T's only once.  perhaps we should repeat it for every starting X?
    X = design(costParams=params['costParams'], fixedParams=params['fixedParams'])
    if params.has_key('Tinit') and params.has_key('Tfinal'):
        Tinit, Tfinal = params['Tinit'], params['Tfinal']
    else:
        Tinit, Tfinal = getTemperatures(X, design, Pinit, Pfinal, params)
    if verbose:
        print 'Tinitial = %.2f, Tfinal = %.2f'%(Tinit,Tfinal)
    if simulated:
        Tfinal = 0.0

    iterationData = []
    for r in xrange(runs):
        if verbose:
            print 'Run #%0.3d,'%r
        sys.stdout.flush()
        X = randomizeInitialX(X)
        runsBestX        = design(X)
        runsBestX.cost   = np.inf

        T = Tinit
        M = Minit
        iter = 0
        while T > Tfinal:
            for i in xrange(int(round(M))):
                try:
                    X = runIteration(X, design, T, params)
                except RuntimeError, e:
                    print e
                    print 'Error: cannot complete run!'
                    print 'Likely b/c network couldn\'t be completed.'
                    X.cost = np.inf
                    break

                if X.cost < runsBestX.cost:
                    runsBestX = design(X)

                iterationData.append([iter, X.cost, runsBestX.cost])
                iter += 1

            T = T/alpha
            M = M*beta

        if verbose:
            print 'Optimal set: ', runsBestX
        sys.stdout.flush()

        if runsBestX.cost < best.cost:
            best = design(runsBestX)
            bestDiscoveries = 1
        elif runsBestX.configParams.values() == best.configParams.values():
            bestDiscoveries += 1
    
    if not verbose2:
        return best
    else:
        return best, iterationData

 
if __name__ == '__main__':
    pass
