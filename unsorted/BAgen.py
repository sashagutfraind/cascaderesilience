'''
Simulations of generation of generalized scale-free networks with variable exponential coefficient

#    Copyright (C) 2007-2010 by  
#    Sasha Gutfraind <ag362@cornell.edu> 
#    Distributed under the terms of the GNU Lesser General Public License 
#    http://www.gnu.org/copyleft/lesser.html 

WARNING:
-----
This implementation is likely not correct

'''

import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
import random
import numpy as np
import numpy.random as npr
import scipy

#import matplotlib
#import matplotlib.pylab as pylab
import pylab
from networkx import *

def barabasi_albert_graph2(n, m, alpha=1, seed=None):
    """
    Variant of the BA model

    A graph of n nodes is grown by attaching new nodes
    each with m edges that are preferentially attached
    to existing nodes with high degree.
   
    :Parameters:
      - `n`: the number of nodes
      - `m`: number of edges to attach from a new node to existing nodes
      - `alpha': probability that attachment is preferential, not random
      - `seed`: seed for random number generator (default=None)

    The initialization is a graph with with m nodes and no edges.

    """
       
    if m < 1 or n < m or alpha < 0 or alpha > 1:
        raise NetworkXError,\
              "NetworkXError must have m>1 and m<n and alpha in [0,1], m=%d,n=%d,alpha=%.3f"%(m,n,alpha)

    if seed is not None:
        random.seed(seed)   

    G=empty_graph(m)       # add m initial nodes (m0 in barabasi-speak)
    G.name="barabasi_albert_graph(%s,%s,%.3f)"%(n,m,alpha)
    edge_targets=range(m)  # possible targets for new edges
    repeated_nodes=[]      # list of existing nodes,
                           # with nodes repeated once for each adjacent edge
    source=m               # next node is m
    while source<n:        # Now add the other n-1 nodes
        numPreferential = npr.binomial(m,alpha)
        G.add_edges_from(zip([source]*numPreferential,edge_targets[1:numPreferential])) # add links to m nodes
        repeated_nodes.extend(edge_targets[1:numPreferential]) # add one node for each new link
        # choose m nodes randomly from existing nodes
        # N.B. during each step of adding a new node the probabilities
        # are fixed, is this correct? or should they be updated.
        # Also, random sampling prevents some parallel edges.

        numRandomEdges = m - numPreferential
        edge_targets2  = random.sample(G.nodes(),numRandomEdges)
        G.add_edges_from(zip([source]*numRandomEdges,edge_targets2)) # add links to m nodes
        repeated_nodes.extend(edge_targets2) # add one node for each new link

        repeated_nodes.extend([source]*m) # and new node "source" has m links
        edge_targets=random.sample(repeated_nodes,m)
        source += 1
    return G


def dorogovtsev_mendes_samukhin(n, m, A=0, seed=None):
    """
    A directed graph based on the directed Dorogovtsev Mendes Samukhin model
        which generalizes the Barabasi-Albert preferential attachment
        to a linear attachment function

    A graph of n nodes is grown by attaching new nodes
    each with m edges that are preferentially attached
    to existing nodes with high degree.
   
    :Parameters:
      - `n`: the number of nodes
      - `m`: number of edges to attach from a new node to existing nodes
      - `A': additional attractiveness, independent of degree
      - `seed`: seed for random number generator (default=None)

    @Article{Dorogovtsev00pref,
      title = {Structure of Growing Networks with Preferential Linking},
      author = {Dorogovtsev, S. N. and Mendes, J. F. F. and Samukhin, A. N.},
      journal = {Phys. Rev. Lett.},
      volume = {85},
      number = {21},
      pages = {4633--4636},
      numpages = {3},
      year = {2000},
      month = {Nov},
      doi = {10.1103/PhysRevLett.85.4633},
      publisher = {American Physical Society},
    }

    """
       
    if m < 1 or n < m or A < 0:
        raise NetworkXError,\
              "NetworkXError must have m>1 and m<n and A in [0,\\infty), m=%d,n=%d,A=%.3f"%(m,n,A)

    if seed is not None:
        random.seed(seed)   

    G=DiGraph()
    G.name="dorogovtsev_mendes_samukhin(%s,%s,%.3f)"%(n,m,A)
    if n < 1:
        return G

    edge_targets = [0] # list of existing nodes that are ends of edges. nodes repeated once for each adjacent edge
                       # initialization with 0 normalizes out
    for newNode in xrange(n):  # Now add the other n-1 nodes
        mx_target    = len(edge_targets)
        new_targets1 = [edge_targets[random.randint(0,mx_target-1)] for x in xrange(mx_target)] 

        G.add_node(newNode)
        nn           = G.number_of_nodes()
        edge_sources = [random.randint(0,nn-1) for x in xrange(m)] 

        new_targets2 = [] #filtered list
        for i,s in enumerate(edge_sources):
            if i == mx_target: 
                break
            t = new_targets1[i]
            if s==t or G.has_edge(s, t):
                continue
            G.add_edge(s,t)
            new_targets2 += [t]
            #iterate, rather than filter since there could be multiple identical new edges
            #TODO: this filtering is useless since in large networks prob. of duplicates is very low
        
        new_targets2.extend(np.floor(A)*[newNode]) 
        if random.random() < A%1.0: 
            new_targets2 += [newNode]
            #TODO: this creates noise.  a partial soln would be to multiply the probability by 10, and also add 10 copies of new_targets2

        edge_targets.extend(new_targets2) 
        #G.remove_edges_from(G.selfloop_edges())
        
    return G

def log_histogram_degrees(G, style):
    degseq=G.degree()
    logdegseq = np.log(degseq)
    
    hist, bins = np.histogram(logdegseq, bins=20)

    L = hist.shape[0]
    slope = slope(np.log(bins), hist, L/3, -2)

    pylab.plot(bins, hists[alpha], style)
    
    return hist, slope

def power_mle(G, xmin=1.):
#based on Clauset et al., http://arxiv.org/abs/0706.1062
    degseq = G.degree()

    #print np.array(degseq).transpose()

    degseqLn = [np.log(xi/xmin) for xi in degseq if xi >= xmin]
    degseqLn.sort() #to reduce underflow. TODO: maybe also normalize and ensure it is summer small to large?

    #print degseqLn
    return 1. + len(degseqLn) / sum(degseqLn)
    
def test_scaling():
    #import matplotlib
    #matplotlib.texmanager.verbose.set_level('debug-annoying')
    #pylab.rc('text', usetex=True)

    #font = {#'family' : 'monospace',
            #'weight' : 'bold',
    #        'size'   : 'smaller',}
    #pylab.rc('font', **font)  # pass in the font dict as kwargs

    slope       = lambda x, y, t0, t1: (y[t1] - y[t0])/(x[t1]- x[t0])
    slopeLogLog = lambda x, y, t0, t1: (np.log(y[t1]) - np.log(y[t0]))/(np.log(x[t1])- np.log(x[t0]))

    nn = 5000
    m  = 10
    #alphas = [0.1, 0.5, 1.0, 2.0, 10.]
    #alphas = [0.2, 1.0, 2.0, 4.0, 20.]
    #alphas = [5., 10., 20.]
    alphas = [1.0]
    graphs = {}
    hists  = {}
    mles   = {} #based on Clauset et al.
    #slopes = {}
    delS   = {}
    styles = dict(zip(alphas, ['r.', 'gs', 'bx']))
    for alpha in alphas:
        #G = barabasi_albert_graph2(nn,m,alpha)
        #G = barabasi_albert_graph(nn,m,alpha)
        #G = dorogovtsev_mendes_samukhin(nn,m,A=alpha).to_undirected()
        G = dorogovtsev_mendes_samukhin(nn,m,A=alpha)
        graphs[alpha] = G

        mles[alpha] = power_mle(G.reverse(), 10.)
        #mles[alpha] = power_mle(G, 10.)
        #hists[alpha], slopes[alpha] = plot_degrees(G)
        
    #pylab.title(r"Performance of the cost evaluations")
    #pylab.xlabel("Number of nodes")
    #pylab.ylabel("Frequency")
    #pylab.figtext(0.15, 0.84, "Power(alpha=0.1: red dots) = %2.3f"%slope[0] + "$\\pm$" +"%2.3f"%delS[0])
    #pylab.figtext(0.2, 0.14, "Power(alpha=0.5: green boxes) = %2.3f"%slope[1] + "$\\pm$" +"%2.3f"%delS[1])
    #pylab.figtext(0.15, 0.80, "Power(alpha=1.0: blue x's) = %2.3f"%slope[2] + "$\\pm$" +"%2.3f"%delS[2])

    #pylab.title(r"\tex[]{Time for a single function evaluation}")

    #print "\nSlopes = %f, %f, %f" %tuple(slopes.values())
    print "\nA\tSlope: " 
    for alpha in alphas:
        print '%.2f:\t%.3f'%(alpha,mles[alpha])

    #pylab.savefig('output/profilecost_results_sparse.eps')
    pylab.show()

    return hists



if __name__ == "__main__":
    test_scaling()
    
