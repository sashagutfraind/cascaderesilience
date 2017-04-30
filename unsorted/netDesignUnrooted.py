'''
Network Generation Programs - unrooted design
    implements the centrality-driven attachment

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
import numpy as np
import numpy.random as npr
import time
import random


#Centrality-driven attachment
class unrootedDesign(netDesign.baseDesign):
    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        if x==None and fixedParams == None:
            self.fixedParams['nn']             = 100 #20
        if x==None and configParams == None:
            self.configParams['edgeFactorA']   = 1.0
            self.configParams['edgeFactorB']   = 0.5  #unless this is raised from 0, no edges might form since ininial cent = 0 and edgeFactorA = 1.0
            self.configParams['maxDegree']     = 3

            self.configParams['pairFactorA']   = 1.0        
            self.configParams['pairFactorB']   = 0.5
            #it might be wiser to do a randomization

    def buildNet(self):
        G                       = nx.Graph()
        nn                      = self.fixedParams['nn']
        G.add_nodes_from(range(nn))

        updatesPerCentrality    = 5
        maxDegree               = self.configParams['maxDegree']
        for update in xrange(nn*5/updatesPerCentrality):
            #print [len(comp) for comp in nx.connected_components(G)]
            centralities = nx.betweenness_centrality(G, normalized=False)
            maxCent      = 1.*max(centralities.values())
            maxCent      = max(maxCent, 1)

            for i in xrange(updatesPerCentrality):
                #if we cannot find a node that has a sufficiently low degree, then we just pick somebody and force a higher degree
                numFailures = 0
                while numFailures < 20:
                    node = npr.randint(0, nn)
                    if len(G.neighbors(node)) < maxDegree:
                        break
                    else:
                        numFailures += 1

                nodeCent = centralities[node]/maxCent
                self.updateNode(G, node, self.configParams, nodeCent)
        
        components = nx.connected_components(G)
        if not len(components) == 1:
            for compNumber in xrange(len(components)-1):
                G.add_edge(random.sample(components[compNumber],   1)[0], 
                             random.sample(components[compNumber+1], 1)[0])
        
        return G

    #def designParamNames(self):
    #    return ['nn','edgeFactorA','edgeFactorB','pairFactorA','pairFactorB','maxDegree']

    def mutate(self, rate=1.0):
        if npr.rand() > rate: return
        #self.hash = None
        self.nets = []
        self.cost = np.inf

        configParams = self.configParams
        param = random.sample(configParams.keys(), 1)[0]
        if param == 'edgeFactorA':
            diff = (npr.rand() - 0.5) * 0.2
            configParams[param] += diff
        elif param == 'edgeFactorB':
            diff = (npr.rand() - 0.5) * 0.2
            configParams[param] += diff
        elif param == 'pairFactorA':
            diff = (npr.rand() - 0.5) * 0.2
            configParams[param] += diff
        elif param == 'pairFactorB':
            diff = (npr.rand() - 0.5) * 0.2
            configParams[param] += diff
        elif param == 'maxDegree':
            diff = round(npr.rand())*2 - 1
            newVal = configParams[param] + diff
            newVal = max(2.0, newVal)
            configParams[param] = newVal

	    return

    def updateNode(self, G, node, configParams, myCent):
        edgeFactorA  = configParams['edgeFactorA']
        edgeFactorB  = configParams['edgeFactorB']
        maxDegree    = configParams['maxDegree']

        pairFactorA  = configParams['pairFactorA']
        pairFactorB  = configParams['pairFactorB']
        
        PrNewEdge = edgeFactorA*myCent + edgeFactorB

        nn = G.number_of_nodes()
        if PrNewEdge > npr.rand() and len(G.neighbors(node)) < maxDegree:
            node2 = npr.randint(0, nn)
            if not G.has_node(node2):
                raise RuntimeError, 'Connecting to a non-existing node!'
            if node != node2:
                G.add_edge(node, node2)
        
        PrPair =    pairFactorA*myCent + pairFactorB

        if PrPair > npr.rand():
            allNbs = G.neighbors(node)

            if len(allNbs) > 1:
                nbs = random.sample(allNbs, 2)
                G.add_edge(nbs[0], nbs[1])

    def vectorEvaluate(self, data_vector):
        #wishlist: test bounds
        self.configParams['edgeFactorA'] = data_vector[0]
        self.configParams['edgeFactorB'] = data_vector[1]
        self.configParams['maxDegree']   = data_vector[2]
        self.configParams['pairFactorA'] = data_vector[3]
        self.configParams['pairFactorB'] = data_vector[4]

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
        data_vector += [self.configParams['edgeFactorA']]
        data_vector += [self.configParams['edgeFactorB']]
        data_vector += [self.configParams['maxDegree']]
        data_vector += [self.configParams['pairFactorA']]
        data_vector += [self.configParams['pairFactorB']]

        return data_vector


if __name__ == '__main__':
    design = unrootedDesign()
    for i in xrange(15):
        design.mutate()
    design.evaluate()
    print design


