'''
Network Generation Programs - G(n,p)

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 
 

'''
import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
#from netDesign import *
import netDesign

import networkx as nx
#from networkx import *
import numpy as np
import numpy.random as npr
import time
import random


#an Erdos-Renyi random graph
class gnpDesign(netDesign.baseDesign):
    #nn = number of nodes
    #p  = connectedness probability between pairs

    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        if x==None and configParams == None:
            self.configParams['p']     = 0.5
            #it might be wiser to do a randomization

    def buildNet(self):
        G                       = nx.Graph()
        nn                      = self.fixedParams['nn']
        G.add_nodes_from(range(nn))

        p = self.configParams['p']
        for node1 in G.nodes():
            for node2 in G.nodes()[node1+1:]:
                if npr.rand() < p:
                    G.add_edge(node1, node2)
        
        return G

    def mutate(self, rate=1.0):
        if npr.rand() > rate: return
        newVal = self.configParams['p'] + (npr.rand() - 0.5) * 0.1
        self.setP(newVal)
        self.purge() #redundant, but just in case

    def setP(self, newVal):
        newVal = max(newVal, 0.0)
        newVal = min(newVal, 1.0)
        self.configParams['p'] = newVal
        self.purge()

    def setParam(self, paramName, newVal):
        if paramName == 'p':
            self.setP(newVal)
        else:
            netDesign.baseDesign.setParam(self, paramName, newVal)

    def vectorEvaluate(self, data_vector):
    #initializes the configuration based on data_vector, evaluates, and returns -c (for use in minimization algorithm)
        #self.configParams['p'] = data_vector[0]might be out of bound
        self.configParams['p'] = 0.5 + 1/np.pi * np.arctan(data_vector[0])
        numTrials = 10
        for trial in xrange(numTrials):
            try:
                self.evaluate()
                break
            except Exception, inst:
                print 'Error has occurred during evaluation.'
                print 'Resubmitting.  trial %d of %d times.'%(trial+1,numTrials)
        return -self.cost
    
    def vectorSave(self):
        data_vector = []
        data_vector.append( np.tan(np.pi*self.configParams['p'] - np.pi/2.) )

        return data_vector

def eff_gnp(n, p, g):
    pass  #can an analytic formula be found for attenuated reach?

def res_gnp(n, p, tau):
    #resilience of gnp, based on Draief et al. JAP 2008.  closely related to res_clique, where tau becomes tau*p
    #WARNING: this code seems to be broken
    '''
    d = netDesignGnp.gnpDesign(fixedParams={'nn':100, 'numNets':10, 'p':0.01})
    d.setParam('tau', 0.1)
    d.evaluate()
    np.average([report['cost']['resilience'] for report in d.netReports])
    netDesignGnp.res_gnp(n=100, p=0.01, tau=0.1)
    '''

    import scipy.optimize
    c     = p*tau*(n-1)
    print c
    if c < 1:
        return 1. - 1./(n-1.) * min(n-1, 1./(1.-c) - 1) #this is only a bound!
    elif c==1.:
        print 'WARNING: No known solution for c=1!'
        if n > 2:
            print 'Computing an upper bound...'
            return 1. - 1./(n-1.) *  (1./(1-p*tau*(n-2)) - 1 )
        else:
            return 1. - 1./(n-1.) *  (1.0)
            #not in the paper, but a sensible estimate
    fnc   = lambda gamma: gamma + np.exp(-c*gamma) - 1. 
    gamma = scipy.optimize.brentq(f=fnc, a=1E-10, b=1.)  #1E-10 is to avoid soln of 0.
    return  1. - gamma


def testConstruction():
    design = gnpDesign()
    for i in xrange(15):
        design.mutate()
    design.evaluate()
    print design
    
if __name__ == '__main__':
    testConstruction()


