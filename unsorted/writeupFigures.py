'''
Various stand-alone figures (for illustrations)

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


'''

import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
import numpy as np
import numpy.random as npr
import random
import matplotlib
#for running scripts silently:
if 'matplotlib.backends' not in sys.modules: matplotlib.use('pdf')
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('PS')
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import networkx as nx
#import pylab
import time, os, sys

import netDesign
import netDesignUnrooted
import netDesignCentrality
import netDesignGnp
import netDesignMulticell

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
#timeNow = lambda : str(time.localtime()[:-6]) + '_' + str(time.localtime()[-6:-3])

def epicycle_layout(G):
    pos = dict.fromkeys(G, None)
    ccs = nx.connected_components(G)
    firstAngle = 0.
    firstRad   = 1.
    secondRad  = 0.2*firstRad
    numCells   = len(ccs)
    for cellNum,cell in enumerate(ccs):
        cellBaseCoords = firstRad*np.array([np.cos(firstAngle),np.sin(firstAngle)])
        secondAngle    = 0.
        for node in cell:
            pos[node] = cellBaseCoords + secondRad*np.array([np.cos(secondAngle),np.sin(secondAngle)])
            secondAngle += 2.*np.pi/len(cell)
        firstAngle += 2.*np.pi/numCells

    return pos

def epicycle_layout_star(G):
    pos = dict.fromkeys(G, None)
    ccs = nx.connected_components(G)
    firstAngle = 0.
    firstRad   = 1.
    secondRad  = 0.2*firstRad
    numCells   = len(ccs)
    for cellNum,cell in enumerate(ccs):
        cellBaseCoords = firstRad*np.array([np.cos(firstAngle),np.sin(firstAngle)])
        secondAngle    = np.pi/3
        leader         = [node for node in cell if G.degree(node) > 1][0]
        pos[leader]    = cellBaseCoords
        cell.remove(leader)
        for node in cell:
            pos[node] = cellBaseCoords + secondRad*np.array([np.cos(secondAngle),np.sin(secondAngle)])
            secondAngle += 2.*np.pi/len(cell)
        firstAngle += 2.*np.pi/numCells

    return pos


def makeFTP_graph(levels=3, add_leaves=True, saveG=False, plot=False):
    random.seed(5)
    npr.seed(5)
    def makeGroup(G, L1node):
        #version with 3 nodes per cell and the leader is extra
        t = '_'+str(random.randint(0,1E5)) 

        G.add_edge(L1node,'L2'+t)
        G.add_edge(L1node,'F11'+t)
        G.add_edge(L1node,'F12'+t)
        G.add_edge(L1node,'F13'+t)
        G.add_edge('F11'+t,'F12'+t)
        G.add_edge('F11'+t,'F13'+t)
        G.add_edge('F12'+t,'F13'+t)

        G.add_edge('L2'+t,'F21'+t)
        G.add_edge('L2'+t,'F22'+t)
        G.add_edge('L2'+t,'F23'+t)
        G.add_edge('F21'+t,'F22'+t)
        G.add_edge('F21'+t,'F23'+t)
        G.add_edge('F22'+t,'F23'+t)

    def makeGroup3c(G, L1node):
        #version with 3 nodes per cell, including the leader
        t = '_'+str(random.randint(0,1E5)) 

        G.add_edge(L1node,'L2'+t)
        G.add_edge(L1node,'F11'+t)
        G.add_edge(L1node,'F12'+t)
        #G.add_edge(L1node,'F13'+t)
        G.add_edge('F11'+t,'F12'+t)
        #G.add_edge('F11'+t,'F13'+t)
        #G.add_edge('F12'+t,'F13'+t)

        G.add_edge('L2'+t,'F21'+t)
        G.add_edge('L2'+t,'F22'+t)
        #G.add_edge('L2'+t,'F23'+t)
        G.add_edge('F21'+t,'F22'+t)
        #G.add_edge('F21'+t,'F23'+t)
        #G.add_edge('F22'+t,'F23'+t)

    G = nx.Graph()
    G.add_edges_from([(0,1),(0,2),(1,2)])
    currentLeaves = [0,1,2]
    newLeaves     = []
    nextNodeNum   = 3
    branching     = 3
    for lev in xrange(levels-1):
        while len(currentLeaves) > 0:
            node          = currentLeaves.pop()
            nodesChildren = []
            for i in xrange(branching):
                G.add_edge(node, nextNodeNum)
                newLeaves    += [nextNodeNum]
                nodesChildren+= [nextNodeNum]
                for syb in nodesChildren:
                    G.add_edge(syb, nextNodeNum)
                nextNodeNum += 1

        currentLeaves = newLeaves
        newLeaves     = []

    
    if add_leaves:
        for leaf in currentLeaves:
            #makeGroup(G, leaf)
            makeGroup3c(G, leaf)

    if saveG:
        #with levels=3: 
        #nx.write_edgelist(G, 'data' + os.sep + 'ftp.graph')
        nx.write_edgelist(G, 'data' + os.sep + 'ftp3c.graph')
        print 'Graph saved'
   
    if plot:
        pylab.figure()
        nx.draw(G)
        pylab.savefig(r'output\ftp.pdf')
    
    return G

def makeCdaFigure():
    d = netDesignCentrality.cdaDesign()
    d.fixedParams['nn'] = 50
    d.configParams['edgeFactorM'] = 9.0
    d.configParams['edgeFactorB'] = 2.0  

    d.fixedParams['transitivityBias'] = 0.9 

    random.seed(5)
    npr.seed(5)
    G = d.buildNet()
    pos = spring_layout2(G, force_ratio=2.)
    pylab.figure()
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    pylab.savefig('output/cda_m=1_b=2_'+timeNow()+'.pdf')

    random.seed(4)
    npr.seed(4)
    d.configParams['edgeFactorM'] = - 3.0
    d.configParams['edgeFactorB'] =  5.0 
    G = d.buildNet()
    pylab.figure()
    pos = spring_layout2(G, force_ratio=2.)
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    pylab.savefig('output/cda_m=-5_b=10_'+timeNow()+'.pdf')

def makeCda2Figure():
    d = netDesignCentrality.cdaDesign2()
    d.fixedParams['nn'] = 50
    d.fixedParams['transitivityBias'] = 0.5 
    d.configParams['edgeFactorM'] = 9.0
    d.configParams['edgeFactorB'] = 2.0  

    random.seed(5)
    npr.seed(5)
    G = d.buildNet()
    pos = spring_layout2(G, force_ratio=2.)
    pylab.figure()
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    #pylab.savefig('output/cda2_'+timeNow()+'.pdf')

    random.seed(5)
    npr.seed(5)
    d.configParams['edgeFactorM'] = - 2.9
    d.configParams['edgeFactorB'] =  8.0 
    G = d.buildNet()
    pylab.figure()
    pos = spring_layout2(G, force_ratio=2.)
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    #pylab.savefig('output/cda2_'+timeNow()+'.pdf')

def makeCliqueFigure():
    d = netDesignMulticell.cliqueDesign()
    d.setParam('nn',42)
    d.setParam('k', 6)
    random.seed(6)
    npr.seed(6)
    d.buildNets()
   
    G = d.nets[0]
    #pos = spring_layout2(G, force_ratio=2.)
    pos = epicycle_layout(G)
    pylab.figure()
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    pylab.savefig('output/clique_'+timeNow()+'.pdf')

def makeCavemenFigure(node_col='red'):
    d = netDesignMulticell.cavemenDesign()
    d.setParam('nn',42)
    d.setParam('k', 6)
    d.setParam('p', 0.5)
    random.seed(9)
    npr.seed(9)  #4, 8
    d.buildNets()
   
    G = d.nets[0]
    pos = spring_layout2(G, force_ratio=1.)
    #pos = epicycle_layout(G)
    pylab.figure()
    nx.draw(G=G, pos=pos, node_color=node_col, with_labels=False, node_size=10)
    pylab.savefig('output/cavemen_'+timeNow()+'.pdf')



def makeCavemenFigure_phaseII(node_col='red'):
    d = netDesignMulticell.cavemenDesign()
    d.setParam('nn',84)
    d.setParam('k', 6)
    d.setParam('p', 0.5)
    random.seed(1)
    npr.seed(1)  #4, 8
    d.buildNets()
   
    G = d.nets[0]
    pos = spring_layout2(G, force_ratio=1.)
    #pos = epicycle_layout(G)
    pylab.figure()
    nx.draw(G=G, pos=pos, node_color=node_col, with_labels=False, node_size=10)
    pylab.savefig('output/cavemen_'+timeNow()+'.eps')


def makeConnectedStarsFigure():
    d = netDesignMulticell.connectedStarsDesign()
    d.setParam('nn',42)
    d.setParam('k', 6)
    d.setParam('p', 0.5)
    #nice seeds: 3, 7
    random.seed(3)
    npr.seed(3)
    d.buildNets()
   
    G = d.nets[0]
    pos = spring_layout2(G, force_ratio=1.)
    #pos = epicycle_layout(G)
    leaders = [node for node in G if G.degree(node) > 1]
    pylab.figure()
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    nx.draw(G=G, pos=pos, node_list=leaders, with_labels=False, node_size=20)
    pylab.savefig('output/connStars_'+timeNow()+'.pdf')
    pylab.savefig('output/connStars_'+timeNow()+'.pdf')

def makeCycleFigure():
    d = netDesignMulticell.cycleDesign()
    d.setParam('nn',42)
    d.setParam('k', 6)
    random.seed(6)
    npr.seed(6)
    d.buildNets()
   
    G = d.nets[0]
    #pos = spring_layout2(G, force_ratio=1.5)
    pos = epicycle_layout(G)
    pylab.figure()
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    pylab.savefig('output/cycle_'+timeNow()+'.pdf')

def makeGnpFigure():
    d = netDesignGnp.gnpDesign()
    d.setParam('p', 0.05)
    d.setParam('nn',42)
    random.seed(4)
    npr.seed(4)
    d.buildNets()
   
    G = d.nets[0]
    pos = nx.spring_layout(G)
    pylab.figure()
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    pylab.savefig('output/gnp_'+timeNow()+'.pdf')

def makeRGG_FigurePhaseII(node_col='blue'):
    G = nx.generators.random_geometric_graph(100, radius=0.15, repel=0.05)
    nx.draw(G=G, pos=G.pos, node_color=node_col, with_labels=False, node_size=10)
    pylab.savefig('output/rgg_'+timeNow()+'.pdf')


def makeStarFigure():
    d = netDesignMulticell.starDesign()
    d.setParam('nn',42)
    d.setParam('k', 6)
    random.seed(6)
    npr.seed(6)
    d.buildNets()
   
    G = d.nets[0]
    #pos = spring_layout2(G, force_ratio=1.5)
    pos = epicycle_layout_star(G)
    leaders = [node for node in G if G.degree(node) > 1]
    pylab.figure()
    nx.draw(G=G, pos=pos, with_labels=False, node_size=10)
    #nx.draw_networkx_nodes(G=G, pos=pos, node_list=leaders, with_labels=False, node_size=10)
    pylab.savefig('output/star_'+timeNow()+'.pdf')


def spring_layout2(G, iterations=50, dim=2, node_pos=None, force_ratio=2.):
    """Spring force model layout
            derived from NetworkX - copyrighted    
    """
    if node_pos==None :  # set the initial positions randomly in 1x1 box
        vpos=nx.random_layout(G, dim=dim)
    else:
        vpos=node_pos
    if iterations==0:
        return vpos
    if G.order()==0:
        k=1.0 * force_ratio
    else:
        k=np.sqrt(1.0/G.order()) * force_ratio# optimal distance between nodes
    disp={}         # displacements

    # initial "temperature" (about .1 of domain area)
    # this is the largest step allowed in the dynamics
    # linearly step down by dt on each iteration so
    # on last iteration it is size dt.
    t=0.1
    dt=0.1/float(iterations+1)
    for i in range(0,iterations):
        for v in G:
            disp[v]=np.zeros(dim)
            for u in G:
                delta=vpos[v]-vpos[u]
                dn=max(np.sqrt(np.dot(delta,delta)),0.01)
                # repulsive force between all
                deltaf=delta*k**2/dn**2
                disp[v]=disp[v]+deltaf
                # attractive force between neighbors
                if G.has_edge(v,u):
                    deltaf=-delta*dn**2/(k*dn)
                    disp[v]=disp[v]+deltaf

        # update positions
        for v in G:
            l=max(np.sqrt(np.dot(disp[v],disp[v])),0.01)
            vpos[v]=vpos[v]+ disp[v]*t/l
        t-=dt
    return vpos

if __name__ == '__main__':
    #make_911_figure()

    #makeCdaFigure()
    #makeCda2Figure()
    #makeCliqueFigure()
    #makeCavemenFigure(node_col='blue')
    #makeConnectedStarsFigure()
    #makeCycleFigure()
    #makeGnpFigure()
    #makeStarFigure()

    #makeCavemenFigure_phaseII(node_col='blue')
    makeRGG_FigurePhaseII(node_col='blue')
    

