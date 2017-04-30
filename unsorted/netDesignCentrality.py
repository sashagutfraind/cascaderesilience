'''
Network Generation Programs - centrality-driven attachment design
    implements the centrality-driven attachment (modified version)

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 

Note:
    This design was rejected because centrality is too sensitive to new edges and too hard to compute
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

#import pdb

#Centrality-driven attachment
#WARNING: this design is not suitable
#1. too fragile when edges are dropped
#2. does not really generate nodes with negative correlation between degree and centrality
class cdaDesign(netDesign.baseDesign):
    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        if x==None and fixedParams == None:
            self.fixedParams['nn']                    = 100
            self.fixedParams['numUpdateCycles']       = 20.
            self.fixedParams['transitivityBias']      = 0.1 
            self.fixedParams['updatesPerCentrality']  = 5. 
        if x==None and configParams == None:
            self.configParams['edgeFactorM']          = 1.0
            self.configParams['edgeFactorB']          = 2.0  
            #unless edgeFactorB is raised from 0, no edges might form 
            #since ininial cent = 0 and edgeFactorA = 1.0
            
            #self.configParams['maxDegree']     = 3
            #self.configParams['pairFactorA']   = 1.0        
            #self.configParams['pairFactorB']   = 0.5
            #it might be wiser to do a randomization

    def buildNet(self):
        G                       = nx.Graph()
        nn                      = self.fixedParams['nn']
        G.add_nodes_from(range(nn))

        #maxDegree               = self.configParams['maxDegree']
        numUpdateCycles         = self.fixedParams['numUpdateCycles']
        updatesPerCentrality    = self.fixedParams['updatesPerCentrality']
        for update in xrange(int(nn*numUpdateCycles/updatesPerCentrality)):
            #print [len(comp) for comp in nx.connected_components(G)]
            centralities = nx.betweenness_centrality(G, normalized=False)
            maxCent      = 1.0*max(centralities.values())
            maxCent      = max(maxCent, 1.0)  #might be < 1

            for i in xrange(int(updatesPerCentrality)):
                node = npr.randint(0, nn)
                ret = self.updateNode(G, node, centralities, maxCent)
                #print ret
        
        return G

    #def designParamNames(self):
    #    return ['nn','edgeFactorA','edgeFactorB','pairFactorA','pairFactorB','maxDegree']

    def mutate(self, rate=1.0):
        if npr.rand() > rate: return
        param = random.sample(self.configParams.keys(), 1)[0]
        if param == 'edgeFactorM':
            diff = (npr.rand() - 0.5) * 0.2
            self.setEdgeFactorM(self.configParams[param] + diff)
        elif param == 'edgeFactorB':
            diff = (npr.rand() - 0.5) * 0.2
            self.setEdgeFactorB(self.configParams[param] + diff)
        self.purge()

    def setEdgeFactorB(self, newVal):
        self.configParams['edgeFactorB'] = newVal
        self.purge()

    def setEdgeFactorM(self, newVal):
        self.configParams['edgeFactorM'] = newVal
        self.purge()

    def setParam(self, paramName, newVal):
        if paramName == 'edgeFactorM':
            self.setEdgeFactorM(newVal)
        elif paramName == 'edgeFactorB':
            self.setEdgeFactorB(newVal)
        else:
            netDesign.baseDesign.setParam(self, paramName, newVal)

    def updateNode(self, G, node, centralities, maxCent):
        #pdb.set_trace()
        myNbs = G.neighbors(node)
        if not self.updateNodeAcceptLink(G.degree(node), centralities[node], maxCent):
            #centrality might rise due to factors beyond this node's control
            if len(myNbs) > 0:
                node2 = random.sample(G.neighbors(node), 1)[0]
                G.delete_edge(node, node2)
                #pdb.set_trace()
                return -1
            else:
                return 0
        
        if npr.rand() < self.fixedParams['transitivityBias'] and len(myNbs) > 0:
            nb_sqred  = set(reduce(lambda x,y: x+y, [G.neighbors(nb) for nb in myNbs]))
            nb_sqred.remove(node)
            if len(nb_sqred) > 0:
                node2     = random.sample(nb_sqred, 1)[0]
            else:
                return 0
        else:
            node2 = npr.randint(0, G.number_of_nodes())
        if not G.has_node(node2):
            raise RuntimeError, 'Connecting to a non-existing node!'
        if node != node2 and self.updateNodeAcceptLink(G.degree(node2), centralities[node2], maxCent):
            G.add_edge(node, node2)
        return 1        

    def updateNodeAcceptLink(self, myDegree, myCent, maxCent):
        edgeFactorM  = self.configParams['edgeFactorM']
        edgeFactorB  = self.configParams['edgeFactorB']
        #maxDegree    = self.configParams['maxDegree']

        #pairFactorA  = self.configParams['pairFactorA']
        #pairFactorB  = self.configParams['pairFactorB']
        
        return myDegree < edgeFactorM*myCent/maxCent + edgeFactorB

    def vectorEvaluate(self, data_vector):
        self.configParams['edgeFactorM'] = data_vector[0]
        self.configParams['edgeFactorB'] = data_vector[1]
        #self.configParams['maxDegree']   = data_vector[2]
        #self.configParams['pairFactorA'] = data_vector[3]
        #self.configParams['pairFactorB'] = data_vector[4]

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
        data_vector += [self.configParams['edgeFactorM']]
        data_vector += [self.configParams['edgeFactorB']]
        #data_vector += [self.configParams['maxDegree']]
        #data_vector += [self.configParams['pairFactorA']]
        #data_vector += [self.configParams['pairFactorB']]

        return data_vector




class centralityTreeDesign(netDesign.baseDesign):
    def __init__(self, x=None, costParams=None, fixedParams=None, configParams=None):
        if x!=None and (costParams!=None or fixedParams!=None or configParams!=None):
            raise ValueError, 'Redundant information in constructor invocation...'

        #self note: this constructor will fail after doing reload() on the base class, unless this class is also reloaded
        netDesign.baseDesign.__init__(self, x=x, costParams=costParams, fixedParams=fixedParams, configParams=configParams)
       
        if x==None and configParams == None:
            self.configParams['b'] = 2.0
            self.configParams['q'] = 1.5  

    def buildNet_forced(self):
        G        = nx.Graph()
        nn       = self.fixedParams['nn']
        currentB = max(self.configParams['b'], 1.)
        q        = self.configParams['q']
        
        G.add_node(0)
        currentLeaves = [0]
        newLeaves     = []
        nextNodeNum   = 1
        while nextNodeNum < nn:
            node = currentLeaves.pop()
            for i in xrange( max(1, int(round(npr.rand()*0.3 + currentB))) ):
                G.add_edge(node, nextNodeNum)
                newLeaves   += [nextNodeNum]
                nextNodeNum += 1
                if nextNodeNum >= nn: 
                    return G

            if currentLeaves == []:
                currentLeaves = newLeaves
                newLeaves     = []
                currentB      = max(currentB*q, 1)  #ensure at least one child is added 
        return G

    def buildNet(self):
        #if the tree cannot be constructed, raises an exception
        G                       = nx.Graph()
        nn                      = self.fixedParams['nn']
        currentB                = self.configParams['b']
        q                       = self.configParams['q']
        
        G.add_node(0)
        currentLeaves = [0]
        newLeaves     = []
        nextNodeNum   = 1
        while nextNodeNum < nn and currentLeaves != []:
            node = currentLeaves.pop()
            for i in xrange(int(round(npr.rand()*0.3 + currentB))):
                G.add_edge(node, nextNodeNum)
                newLeaves   += [nextNodeNum]
                nextNodeNum += 1
                if nextNodeNum >= nn: 
                    return G

            if currentLeaves == []:
                currentLeaves = newLeaves
                newLeaves     = []
                currentB      = currentB*q 
                #ensure at least one child is added 
                #currentB      = max(currentB*q, 1)
        if nextNodeNum != nn:
            raise ValueError, 'Tree cannot be constructed with %d nodes.  Increase b or decrease q'%nn
            #for node in xrange(nextNodeNum, nn+1):
            #    G.add_node(node)
        return G


    def mutate(self, rate=1.0):
        if npr.rand() > rate: return
        param = random.sample(self.configParams.keys(), 1)[0]
        if param == 'b':
            diff = (npr.rand() - 0.5) * 0.2
            self.setB(self.configParams[param] + diff)
        elif param == 'q':
            diff = (npr.rand() - 0.5) * 0.2
            self.setQ(self.configParams[param] + diff)
        self.purge()

    def setParam(self, paramName, newVal):
        if paramName == 'b':
            self.setB(newVal)
        elif paramName == 'q':
            self.setQ(newVal)
        else:
            netDesign.baseDesign.setParam(self, paramName, newVal)

    def setB(self, newVal):
        self.configParams['b'] = max(newVal, 1.)
        self.purge()

    def setQ(self, newVal):
        self.configParams['q'] = newVal
        self.purge()

    def vectorEvaluate(self, data_vector):
        self.setParam('b', data_vector[0])
        self.setParam('q', data_vector[1])

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
        data_vector += [self.configParams['b']]
        data_vector += [self.configParams['q']]

        return data_vector



def centrality_correlation(G=None):
    import scipy, pylab, cPickle
    if G==None:
        npr.seed(3)
        random.seed(3)
        des = centralityTreeDesign()
        #des = cdaDesign()

        #gives high correlation
        des.setParam('b', 15.)
        des.setParam('q', 0.7)

        #gives low correlation
        #des.setParam('b', 2.0)
        #des.setParam('q', 3.0)
        des.fixedParams['nn'] = 300
        G   = des.buildNet()
    
    results      = []
    centralities = nx.centrality.brandes_betweenness_centrality(G)
    for node in G:
        if G.degree(node) == 1: 
            continue
        #jitter for visualization:
        #results.append((G.degree(node)+npr.rand()*0.4, centralities[node]))
        
        results.append((G.degree(node), centralities[node]))

    results = np.array(results)
    pylab.rc('text', usetex=True)
    pylab.plot(results[:,0], results[:,1], '.')
    #pylab.title(r'lambda=%2.1f'%(rat,))
    pylab.xlabel('Degree')
    pylab.ylabel('Betweenness')
    
    #fixme: statistics question - what are all the other coefficients:?
    corr = scipy.corrcoef(results[:,0], results[:,1])[0,1]
    pylab.figtext(0.2, 0.8, 'Correlation=%f'%corr)
    
    filename   = 'output/correlation_motion.vs.betweenness.pkl'
    outputFile = open(filename, 'wb')
    report = {'results':results, 'correlation':corr}
    cPickle.dump(report, outputFile)
    outputFile.close()
    print 
    print 'Pickle: ' + filename + ' written!'

    results = np.array(results)
    print 'Correlation: %.4f'%corr

    try:
        pylab.savefig('output/correlation_motion.vs.betweenness_lambda='+str(rat)+'_results.eps')
    except:
        print 'Unable to save figure...'



if __name__ == '__main__':
    design = unrootedDesign()
    for i in xrange(15):
        design.mutate()
    design.evaluate()
    print design


