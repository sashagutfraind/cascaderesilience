'''
computes the efficiency of a graph (various metrics)

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


Notes:
-----

'''

import sys, pdb
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
import numpy as np
import numpy.random as npr
import networkx as nx
#from networkx import *

def avgDistance(G):
    #1. doesn't work on disconnected graphs
    #2. inefficient: the computation should actually take C(n,2) rather than n^2 in time
    return np.average(nx.eccentricity(G, with_labels=True).values())

def invDiameter(G):
    if nx.is_connected(G):
        return 1.0/nx.diameter(G)
    else:
        return 0.0

def ringReach(G, params={}, normalize=False):
    '''computes the average size of a node's distance to other nodes, attenuated by distance to power g
    preconditions:
        nodes are indexed 0..nn-1
    '''
    sys.path.append('/home/gfriend/lib/python')
    sys.path.append('/homes/gfriend/lib/python')
    g = params.get('g', 1.0)
    weighted = params.get('weighted', False)

    nn = G.number_of_nodes()
    if nn == 0:
        print 'Warning: ring reach undefined on empty graph! Returning 0.0'
        return 0.0

    counts = np.zeros(nn, dtype=np.double) #help deal with underflow

    if not weighted:
        distFunc = nx.single_source_shortest_path_length  #BFS
        wtNorm = 1.
    else:
        distFunc = nx.single_source_dijkstra_path_length
        #scipy.stats.hmean
        invWts = [(1./(edge[2]['weight']))**g for edge in G.edges(data=True)]
        invWts.sort()  #underflow control
        wtNorm = len(invWts) / sum(invWts)

    for node in G:
        nodeSum = np.double(0.)
        dist_to_n = distFunc(G,node)
        dist_to_n.pop(node)
        for d in dist_to_n.values():
            nodeSum += (1./d)**g
        counts[node] = nodeSum

    counts.sort() #help with underflow 
    ret = counts.sum()

    #this normalization is problematic for weighted case: pair distances could be shorter than 1.
    if normalize:
        return ret * wtNorm / (nn*(nn-1))
    else:
        return ret

def testRingReach():
    def compHelper(val, expVal, tolerance=0.01):
        print 'Computed: %f'%val
        print 'Expected: %f'%expVal
        print
        assert abs(val-expVal) < tolerance
        
    print 'Testing Ring Reach:'
    G = nx.complete_graph(10)
    exp = 10 * 9
    compHelper(ringReach(G=G, params={'g':2.0}, normalize=False), exp)

    G = nx.star_graph(10)  #11 nodes
    exp = 1*10.0 + 10*1.0  #+small change
    compHelper(ringReach(G=G, params={'g':100.0}, normalize=False), exp)

    G = nx.Graph()
    G.add_edge(0,1,weight=1)
    G.add_edge(0,2,weight=2)
    G.add_edge(1,2,weight=2)
    compHelper(ringReach(G=G, params={'g':1.0, 'weighted':True}, normalize=True), 1.0)

    G = nx.Graph()
    G.add_edge(0,1,weight=1)
    G.add_edge(0,2,weight=2)
    compHelper(ringReach(G=G, params={'g':1.0, 'weighted':True}, normalize=True), 44./54)

    print 'Test passed!'

if __name__ == '__main__':
    testRingReach()

