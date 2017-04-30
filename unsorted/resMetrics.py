'''
Computes the resilience of a graph to an epidemic

    Copyright (C) 2007-2013 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


Notes:
-----

includes:
    various implementations of various measures;  invokes resMetricsX.py
    tests

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

def caseCounterGen():
    i = 0
    while True:
        yield str(i)
        i += 1

def demoKeeling():
    import pylab
    tau = 0.5
    params = {'tau':tau}

    G = nx.generators.newman_watts_strogatz_graph(100,4,0.13,1)
    #G = generators.erdos_renyi_graph(100,0.1,1)
    #G = generators.cycle_graph(100)
    #times, recovered = resilienceKeeling(G, params={}) 
    #pylab.plot(times, recovered, 'r.')
    times, infected = resilienceKeeling(G, params={'tau':9.1, 'g':0.1})  #how can tau>1.0?
    pylab.plot(times, infected, 'b.')
    pylab.show()

   #todo time


def extentKeelingDx(X,t,tau,g,phi,nu,nn):
#the vector field for the solver
    def abc(phi,nu,nn,ab,bc,ac,a,b,c): 
        if a < 1E-30 or b < 1E-30 or c < 1E-30:
           return 0.
        else:
           return 1.*(nu-1)/nu*ab*bc/b*( 1-phi + phi*nn/nu*ac/a/c )

    ss  = X[0]
    si  = X[1]
    sr  = X[2]
    ii  = X[3]
    ir  = X[4]

    s   = (2*ss+  si+sr)/nu
    i   = (  si+2*ii+ir)/nu
    r   = nn-s-i
    #print s,i,r

    ssi = abc(phi,nu,nn,ss,si,si,s,s,i)
    isi = abc(phi,nu,nn,si,si,ii,i,s,i)
    rsi = abc(phi,nu,nn,sr,si,ir,r,s,i)
    ssdot = -2*tau*ssi
    sidot = tau*(ssi-isi-si) - g*si
    srdot = -tau*rsi + g*si
    iidot = 2*tau*(isi+si) - 2*g*ii
    irdot = tau*rsi+g*(ii-ir)

    return (ssdot,sidot,srdot,iidot,irdot)


def extentKeelingGraphPhi(G):
#computes the phi value of the graph (for Keeling)
    mat   = nx.to_numpy_matrix(G)
    matsq = mat*mat
    matq  = matsq*mat

    triangles = matq.trace()[0,0]
    triplets  = matsq.sum() - matsq.trace()[0,0]
    return 1.*triangles/triplets


def extentKeeling(G, params={}):
    #based on Keeling 1999: Effects of local spatial structure on epidemiological invasions
    import scipy.integrate

    tau = params.get('tau', 0.1) 
    g   = params.get('g', 0.1)
    phi = extentKeelingGraphPhi(G)
    nu  = np.average(degree(G))
    nn  = G.number_of_nodes()*1.
    ne  = G.number_of_edges()*1.

    times = np.array([0, 1./g*ne*100])

    ss  = -1.*nu + 1.*ne
    si  = +1.*nu
    sr  = 0.
    ii  = 0.
    ir  = 0.

    r   = 0.
    params = (tau,g,phi,nu,nn)
    '''
    for trial in range(10):
        x0     = (ss,si,sr,ii,ir)
        solns  = scipy.integrate.odeint(keelingDx,x0,times,params)
        ss  = solns[-1,0]
        si  = solns[-1,1]
        sr  = solns[-1,2]
        ii  = solns[-1,3]
        ir  = solns[-1,4]
        rprime = nn - (2*ss + 2*si + sr + 2*ii + ir)/nu
        if (r-rpime)/nn < 0.001:
           r = rprime
           break
        r = rprime

    return r
    '''
    #times = np.arange(0,1./g*nn*10,g/10)
    times = np.arange(0,tau/g,g/30)
    x0    = (ss,si,sr,ii,ir)
    solns = scipy.integrate.odeint(extentKeelingDx,x0,times,params)
    #return times, np.repeat(np.double(ne),len(times))-solns[:,0]-solns[:,1]-solns[:,3]
    
    print 'Graph: avg degree=%2.2f, triangles/triplets=%2.2f'%(nu,phi)
    #wishlist: report expected epidemic threshold and the computed value
    #wishlist: verify that this returns the absolute number of infected nodes, not fraction
    return times, (2*solns[:,1]+solns[:,3]+solns[:,4])/nu


def extentSIRSimL(G, params={}):
    #list-dict implementation of the contagion. much faster than a queue or a list-list
    #lavishly expensive: 5 repetitions per G.node
    sys.path.append('/home/gfriend/lib/python')
    sys.path.append('/homes/gfriend/lib/python')
    nn          = G.number_of_nodes()
    tau		    = params.get('tau', 0.5)
    repetitions = params.get('repetitions', 5)
    if 'tau' not in params:
        print 'Warning: using detault tau!'

    if G.nodes() == []: return 0., 0.

    extents = np.zeros(nn*repetitions, dtype=np.double)
    for startN in G.nodes():
        for r in xrange(repetitions):
            extents[r+startN*r] += extentSIRSim_helper(G=G, startN=startN, params={'tau':tau}, vaccinated=[], weighted=False)
            
    return np.average(extents), np.std(extents)


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


############################

extentSIRSim = extentSIRSimLAdaptive

############################


def extentSIRclique(G, params={}):
    nn          = G.number_of_nodes()
    tau 	    = params.get('tau', 0.5)
    g           = params.get('tau', 0.1)
    repetitions = params.get('repetitions', 5)

    pass


def extentSIRtreeAnalytic(G, params={}):
#note: this should be extended to consider value of nodes, if possible

    #wishlist: compute the variance
    tau = params.get('tau', 0.5)

    degrees      = nx.degree_histogram(G)
    ladder       = np.arange(0, len(degrees))
    normP        = 1.0*G.number_of_nodes()
    G_null_prime = np.dot(ladder, degrees) / normP

    avgDeg       = G_null_prime
    ladder2      = np.concatenate(([0], ladder[:-1]))
    G_one_prime  = np.dot(ladder2, degrees*ladder) / normP / avgDeg 
    
    #formula (22) in Newman, 2002
    mean  = 1 + (tau*G_null_prime)/(1-tau*G_one_prime)
    
    #fixme
    stdev = np.NaN 

    #wishlist: verify that this returns the absolute number of infected nodes, not fraction
    return mean, stdev


def profileCythonSim():
    params = {'tau':0.1, 'tolerance':0.05, 'minReplications':100}

    print 'Python implementation:'
    G = nx.generators.erdos_renyi_graph(100,0.3)
    t = timeit.default_timer() 
    print extentSIRSimLAdaptive(G, params=params)
    print timeit.default_timer()-t

    print 'Cython implementation:'
    t = timeit.default_timer() 
    try:
      print extentSIRSimX(G, params=params)
    except Error, e:
      print e
    print timeit.default_timer()-t


def testAnalytic():
    def reportHelper(G, params, expected=None):
        meanNum = extentSIRSim(G, params)[0]
        meanCyt = extentSIRSimX(G, params=params)[0]

        print 'Mean Python: %f'%meanNum
        if expected!=None and abs(meanNum-expected) > 1:
            print 'Test failed: displacement > 1 node'
        else:
            print 'Test passed!'
        print 'Mean Python: %f'%meanCyt
        if expected!=None and abs(meanCyt-expected) > 1:
            print 'Test failed: displacement > 1 node'
        else:
            print 'Test passed!'
            

    tau = 0.5
    params = {'tau':tau, 'tolerance':0.01, 'minReplications':100}

    print 'Params: '+str(params)
    print

    BT = nx.generators.balanced_tree(10,2)
    print 'Balanced tree 10,2: ' 
    reportHelper(G=BT, params=params, expected=None)
    
    C10 = nx.generators.cycle_graph(10)
    print 'Cycle 10: '
    reportHelper(G=C10, params=params, expected=None)
    
    C100 = nx.generators.cycle_graph(100)
    print 'Cycle 100: '
    reportHelper(G=C100, params=params, expected=None)

    #expected for a large cycle:
    #1*(1-tau)**2 + 2*2*tau*(1-tau)**2 + 3*3*tau*tau*(1-tau)**2 + 4*4*tau*tau*tau*(1-tau)**2 + 
    s = 0.
    for i in xrange(1,100):
        s += (i**2) * (tau**(i-1)) * ((1-tau)**2)
    print
    print 'For a cycle, expect: %f'%s

    print extentSIRSimL(C10, params)
    #print extentSIRSimL(G, params)
    #print extentSIRSimQ(G, params)


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
            m, s  = extentSIRSim(G, params)
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



def testSequetialSampling():
    def exBruteSampler(G, params):
        nn = G.number_of_nodes()
        numTrials = nn*params['repetitions']
        
        sum = 0.
        print 'Num trials: %d'%numTrials
        for i in xrange(numTrials):
            sum += npr.rand()

        return sum / float(numTrials)

    def exSequentialSampler(G, params):
        #nn = G.number_of_nodes()
        minReplications = params.get('minReplications', 20)
        tolerance       = params.get('tolerance', 0.1)
        
        outcomes   = []
        #np.std hould use ddof=1, but it is not compatible with all versions
        acceptable = lambda outcomes, tolerance: 1.96*np.std(outcomes) / np.sqrt(len(outcomes)) < tolerance
        trial      = 1
        while trial <= minReplications or not acceptable(outcomes, tolerance):
            outcomes += [npr.rand()]
            trial    += 1
            #print outcomes[-1]

        print 'Num trials: %d'%(trial+1)
        return np.average(outcomes)

    G = nx.generators.erdos_renyi_graph(100, 0.1)
    
    paramsBrute={'tau': 0.1, 'repetitions':5}
    print 'Brute force:'
    print extentSIRSimL(G, params=paramsBrute)
    #print exBruteSampler(G, params=paramsBrute)
    
    paramsAdaptive=paramsBrute.copy()
    paramsAdaptive['minReplications'] = 40
    paramsAdaptive['tolerance']       = 0.01
    #minReplications should be kept to at least 10
    #tolerance of 0.1 seems to be too low, but 0.01 is too high; 
    #note: n grows quadratically in tolerance

    print 'Adaptive:'
    print extentSIRSimLAdaptive(G, params=paramsAdaptive)
    #print exSequentialSampler(G, params=paramsAdaptive)
    



if __name__ == '__main__':
    print 'Running Tests...\n'
    profileCythonSim()
    testCython()
    #testAnalytic()
