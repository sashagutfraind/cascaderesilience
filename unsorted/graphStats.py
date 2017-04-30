'''
Calculation of various graph properties

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


'''

import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')

import networkx as nx
import numpy as np
import numpy.random as npr
import random


def avgDegree(G):
    return np.average(nx.degree(G).values())

def medDegree(G):
    return np.median(nx.degree(G).values())

def maxDegree(G):
    return np.max(nx.degree(G).values())

'''
def avgDegree(G):
#better just np.average(nx.degree(G))
    degrees      = nx.degree_histogram(G)
    ladder       = np.arange(0, len(degrees))
    normP        = 1.0*G.number_of_nodes()
    return np.dot(ladder, degrees) / normP
'''

def sMetric(G, normalized=True):
#degree correlations, as in Willinger et al.
    D   = np.matrix(nx.degree(G))
    mat = nx.to_numpy_matrix(G)
    
    val = 0.5 * ((D*mat*D.T)[0,0]) 
    if normalized:
        val /= (np.multiply(D,D)*D.T)[0,0] * 0.5

    return val

    
if __name__ == '__main__':
    print 'Done'

