'''
Computes the resilience of a graph to an epidemic (implementation in Cython)

    Copyright (C) 2007-2013 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


Notes:
-----

'''

import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
#sys.path.append('./extensions')

import numpy as np
cimport numpy as np
import numpy.random as npr

cdef extern from "stdlib.h":
    #TODO: random,srandom fails to compile under mingw.  does new code work in linux?
    #long c_libc_random "random"()
    #void c_libc_srandom "srandom"(unsigned int seed)
    long c_libc_random "rand"()
    void c_libc_srandom "srand"(unsigned int seed)
        
from stdlib cimport *

DEF MAX_SEED = 32767


############################

extentSIRSim = extentSIR_adaptive

############################


cpdef extentSIR_adaptive(G, params={}):
#main method: a wrapper
    cdef int nn     = G.number_of_nodes()
    cdef float tau	= params.get('tau', 0.5)
    cdef int minReplications = params.get('minReplications', max(nn/10, 40))
    cdef float tolerance = params.get('tolerance', 0.05)
    if 'minReplications' not in params:
        print 'Warning: using detault minReplications=%d!'%minReplications
    if 'tau' not in params:
        print 'Warning: using detault tau=%.3f!'%tau
    if 'tolerance' not in params:
        print 'Warning: using detault tolerance=%.3f!'%tolerance

    if nn==0: return 0., 0.  #this is a problem since we assume initial infection site
    elif nn==1: return 1., 0. 

    c_libc_srandom(npr.randint(low=MAX_SEED)) #could potentially pass the seed and avoid importing npr
    
    cdef int **lol = <int **>malloc(nn * sizeof(int*)) #graph as list-of-lists
    cdef int node,
    cdef list nbs
    cdef int numNbs
    cdef int j
    for node in range(nn): #range is apparently better than xrange
        nbs = G.neighbors(node)
        #print nbs
        numNbs = len(nbs)
        lol[node] = <int *>malloc((1+numNbs) * sizeof(int))
        lol[node][0] = numNbs #0 stores number of neighbors
        for j from 1<=j<=numNbs:
            lol[node][j] = nbs[j-1]


    cdef float* ret= _extentSIR_adaptive(lol, nn, tau, minReplications, tolerance)
    for node in range(nn):
        free(lol[node])
    free(lol)
    cdef float mean = ret[0]
    cdef float std  = ret[1]
    if mean > nn: 
        raise RuntimeError, 'Something went wrong: mean size > nn'
    if mean < 1.: 
        raise RuntimeError, 'Something went wrong: mean size < 1.'
    return mean, std

cdef float* _extentSIR_adaptive(int** lol, int nn, float tau, int minReplications, float tolerance):
#helper method in Cython
    cdef list extents = []  #using lists is slow but switching to ndarray is not improving things
    cdef int trial   = 1
    cdef int startNode
    while (trial <= minReplications or not stdAcceptable(extents, float(nn), tolerance)): # and trial < maxReplications:
        startNode = c_libc_random() % nn
        #print
        #print startNode
        extents  += [runSingleNode(lol, startNode, tau, nn)]
        trial    += 1

    #print 'Num trials: %d'%trial 
    cdef float ret[2]
    ret[0] = np.average(extents)
    ret[1] = np.std(extents)

    return ret

cdef int runSingleNode(int** lol, int startN, float tau, int nn):
    cdef int intThreshold = <int>(MAX_SEED*tau)
    cdef np.ndarray[np.int_t, ndim=1] uninfectable = np.zeros(nn, dtype=int)
    uninfectable[startN] = 1
   
    cdef np.ndarray[np.int_t, ndim=1] infected    = np.zeros(nn, dtype=int)
    infected[0] = startN
    cdef int head = 0
    cdef int tail = 1
    cdef int nb
    cdef int nbCounter
    cdef int numNbs
    while head < tail:
        curNode = infected[head]
        #print uninfectable
        #print infected
        #print curNode, head, tail
        head   += 1
        numNbs = lol[curNode][0]
        for nbCounter from 1<=nbCounter<=numNbs:
            nb = lol[curNode][nbCounter]
            if uninfectable[nb]==0 and intThreshold > (c_libc_random()%MAX_SEED):
                #print 'new inf'
                uninfectable[nb] = 1
                infected[tail]   = nb
                tail            += 1
            #else:
            #    print 'inf failed'
            #    if uninfectable[nb]!=0: print 'nb already infected'
    return tail

cdef int stdAcceptable(list outcomes, float norm, float tolerance):
    cdef float val = 1.96*np.std(outcomes) / norm / np.sqrt(len(outcomes))
    #print val
    return val < tolerance


