'''
Runs various epidemic models

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
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('PS')
if 'matplotlib.backends' not in sys.modules: matplotlib.use('pdf')
from matplotlib import rc
rc('text', usetex=True)

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import pylab
import networkx as nx  #must be after pylab, else gives problems with X
import time, os, sys
import csv
import cPickle
import pp
import taskDispatcher
import pdb

#import graphStats
#import netDesign
import resMetrics

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
#timeNow = lambda : str(time.localtime()[:-6]) + '_' + str(time.localtime()[-6:-3])

#    Copyright (C) 2004 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

import string

def compare_stochastic_threshold():
    clique = nx.complete_graph(40)
    clique.name = 'K%d'%clique.number_of_nodes()
    G = clique
    
    #grid = nx.convert_node_labels_to_integers(nx.grid_2d_graph(10,5))
    #grid.name = 'grid10x5'

    #trirand = triangle_rand_g()
    #graphs = [         diamond, ]
    #lams = [1, 1.5, 2.0, 2.5]
    #lams = [1., 2., 3., 4.]
    #lams = [2., 4., 6., 8.]
    lams = [1., 2., 4., 8., 12., 16.]
    #lams = [10., 20., 30., 40.]
    Ks    = []
    for lam in lams:
        Kvals = dict([(node,npr.poisson(lam)) for node in G])
        Ks.append(Kvals)
    names = ['$\\lambda=%.1f$'%lam for lam in lams]
    markers = ['r-^','b-o', 'c-', 'g-', 'y-.', 'r.-', 'k-.', 'b-.']
    linewidths = [1,2,3,4,5,6,7,8]
    pylab.figure()
    vsTau = lambda Kvals, tau: resMetrics.extentSIRSimLAdaptive(G, 
                                             params={'tau':tau, 'vaccinated':[], 
                                             'transmissionModel':'KFleeSIR', 'weighted':False, 
                                             'minReplications':40, 'tolerance':0.005, 'K':Kvals})[0]
    pylab.hold(True)
    taus = np.arange(0, 1.01, 0.05)
    for i,Kvals in enumerate(Ks):
        extentFunc = lambda tau: vsTau(tau=tau, Kvals=Kvals)
        extents = np.array(map(extentFunc, taus))
        pylab.plot(taus,extents,markers[i],linewidth=linewidths[i],label=G.name+'(%s)'%names[i])
    pylab.hold(False)
    #pylab.legend(loc='best')
    pylab.legend(loc='upper center')
    pylab.xlabel(r'$\alpha$', fontsize=20)
    #pylab.ylabel('fraction of nodes infected', fontsize=20)
    pylab.ylabel('Number of nodes infected', fontsize=20)
    pylab.savefig('output/epidemics_poissonK_'+timeNow())
    pylab.show()



def configuration_model():
    func = lambda n:nx.utils.powerlaw_sequence(n, 4.0)
    z=nx.create_degree_sequence(100, func)
    G=nx.configuration_model(z)
    G=nx.Graph(G)
    G.remove_edges_from(G.selfloop_edges())
    G.name = 'power law'
    
    return G
    
def davis_club_graph(create_using=None, **kwds):
    nwomen=14
    nclubs=18
    G=nx.generators.empty_graph(nwomen+nclubs,create_using=create_using,**kwds)
    G.clear()
    G.name="Davis Southern Club Women"

    women="""\
EVELYN
LAURA
THERESA
BRENDA
CHARLOTTE
FRANCES
ELEANOR
PEARL
RUTH
VERNE
MYRNA
KATHERINE
SYLVIA
NORA
HELEN
DOROTHY
OLIVIA
FLORA"""

    clubs="""\
E1
E2
E3
E4
E5
E6
E7
E8
E9
E10
E11
E12
E13
E14"""

    davisdat="""\
1 1 1 1 1 1 0 1 1 0 0 0 0 0
1 1 1 0 1 1 1 1 0 0 0 0 0 0
0 1 1 1 1 1 1 1 1 0 0 0 0 0
1 0 1 1 1 1 1 1 0 0 0 0 0 0
0 0 1 1 1 0 1 0 0 0 0 0 0 0
0 0 1 0 1 1 0 1 0 0 0 0 0 0
0 0 0 0 1 1 1 1 0 0 0 0 0 0
0 0 0 0 0 1 0 1 1 0 0 0 0 0
0 0 0 0 1 0 1 1 1 0 0 0 0 0
0 0 0 0 0 0 1 1 1 0 0 1 0 0
0 0 0 0 0 0 0 1 1 1 0 1 0 0
0 0 0 0 0 0 0 1 1 1 0 1 1 1
0 0 0 0 0 0 1 1 1 1 0 1 1 1
0 0 0 0 0 1 1 0 1 1 1 1 1 1
0 0 0 0 0 0 1 1 0 1 1 1 1 1
0 0 0 0 0 0 0 1 1 1 0 1 0 0
0 0 0 0 0 0 0 0 1 0 1 0 0 0
0 0 0 0 0 0 0 0 1 0 1 0 0 0"""


    # women names
    w={}
    n=0
    for name in women.split('\n'):
        w[n]=name
        n+=1

    # club names
    c={}
    n=0
    for name in clubs.split('\n'):
        c[n]=name
        n+=1

    # parse matrix
    row=0
    for line in davisdat.split('\n'):
        thisrow=list(map(int,line.split(' ')))
        for col in range(0,len(thisrow)):
            if thisrow[col]==1:
                G.add_edge(w[row],c[col])
        row+=1
    # return graph and women and clubs lists
    (G,women,clubs) = (G,list(w.values()),list(c.values()))

    def project(B,pv,result=False,**kwds):
        """
        Returns a graph that is the unipartite projection of the
        bipartite graph B onto the set of nodes given in list pv.

        The nodes retain their names and are connected if they share a
        common node in the vertex set of {B not pv}.

        No attempt is made to verify that the input graph B is bipartite.
        """
        if result:
            G=result
        else:
            G=nx.Graph(**kwds)
        for v in pv:
            G.add_node(v)
            for cv in B.neighbors(v):
                G.add_edges_from([(v,u) for u in B.neighbors(cv)])
        return G

    # project bipartite graph onto women nodes
    W=project(G,women)
    W = nx.convert_node_labels_to_integers(W)
    W.name = 'southern women'

    return W

def karate_graph(create_using=None, **kwds):
    from networkx.generators.classic import empty_graph

    G=empty_graph(34,create_using=create_using,**kwds)
    G.name="Karate Club"

    zacharydat="""\
0 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0
1 0 1 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 0 1 0 0 0
1 1 0 1 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 0
1 1 1 0 0 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 0 1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 1
0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 1
0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 1 0 0 1 0 1 0 1 1 0 0 0 0 0 1 1 1 0 1
0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 1 0 0 1 1 1 0 1 1 0 0 1 1 1 1 1 1 1 0"""


    row=0
    for line in zacharydat.split('\n'):
        thisrow=list(map(int,line.split(' ')))
        for col in range(0,len(thisrow)):
            if thisrow[col]==1:
                G.add_edge(row,col) # col goes from 0,33
        row+=1
    return G




def plot_real():
    karate=karate_graph()
    women =davis_club_graph()
    #graphs = [women, karate,]
    graphs = [       karate,]
    names = ['SIR', '2FleeSIR']
    markers = ['r-^','b-o']
    linewidths = [1,2,3,4]
    pylab.figure()
    vsTauG = lambda G, tau, transmissionModel: resMetrics.extentSIRSimLAdaptive(G, params={'tau':tau, 'vaccinated':[], 'transmissionModel':transmissionModel, 'weighted':False, 'minReplications':40, 'tolerance':0.003})[0]
    for j,G in enumerate(graphs):
        funcs = [lambda tau: vsTauG(G, tau, 'SIR'), lambda tau: vsTauG(G, tau, 'KFleeSIR')]
        taus = np.arange(0, 1.01, 0.05)
        pylab.hold(True)
        for i in xrange(2):
            extents = np.array(map(funcs[i], taus))
            pylab.plot(taus,extents,markers[i],linewidth=linewidths[i],label=G.name+'(%s)'%names[i])
    pylab.hold(False)
    pylab.legend(loc='best')
    pylab.xlabel(r'$\tau$', fontsize=20)
    pylab.ylabel('Number of nodes infected', fontsize=20)
    pylab.savefig('output/epidemics'+timeNow())
    pylab.show()

def plot_sirf(fname=None):
    alphaEli = []
    extentEli = [] 
    with open('data/eli_noFleeStoch.dat', 'rb') as f:
        dataDict = csv.DictReader(f, delimiter=' ')
        for record in dataDict:
            alphaEli.append(float(record['alpha']))
            extentEli.append(float(record['fraction_final_infected']))
    
    reports = {}
    if fname == None:
        #fname = 'sirf/meanExtent_NumNodes=1000_2011_3_15__4-2-1628622592_.csv'
        #fname = 'sirf/meanExtent_NumNodes=10000__2011_3_12__21-8-4198687.csv'
        #fname = 'sirf/meanExtent_NumNodes=1000__2011_3_12__35-4-4198447.csv'
        fname = 'sirf/meanExtent_NumNodes=10000__2011_4_17__15-7-12.csv'
        #fname = 'sirf/meanExtent_NumNodes=1000000__2011_4_17__17-37-11.csv'
        numNodes = 10000.
    with open(fname, 'rb') as f:
        dataDict = csv.DictReader(f, delimiter=',')
        for record in dataDict:
            alpha   = float(record['alpha'])
            gamma   = float(record['gamma'])
            extent  = float(record['meanExtent'])/numNodes
            reports[(alpha,gamma)] = extent
    alphas = list(set([alpha for alpha,gamma in reports]))
    alphas.sort()
    #alphas = alphas[:20]  #FIXME
    gammas = list(set([gamma for alpha,gamma in reports]))
    gammas.sort()
    #gammas = gammas[:20]  #FIXME

    X, Y = np.meshgrid(gammas, alphas)
    #Z    = np.apply_along_axis(lambda (x,y): reports[(x,y)], 0, (X,Y))
    Z    = np.zeros(X.shape)
    for i, alpha in enumerate(alphas):
        for j, gamma in enumerate(gammas):
            Z[i,j] = reports[(alpha,gamma)]
    
    fig=plt.figure()
    plt.hold(True)
    plt.plot(alphaEli, extentEli, label='gamma=0 (eli)')
    plt.plot(alphas, Z[:,gammas.index(0.)], label='gamma=0')
    plt.plot(alphas, Z[:,gammas.index(100.)], label='gamma=100')
    plt.plot(alphas, Z[:,gammas.index(1000.)], label='gamma=1000')
    plt.plot(alphas, Z[:,gammas.index(10000.)], label='gamma=10000')
    #print alphas
    #print Z[:,0]
    plt.xlabel('alpha')
    plt.ylabel('attack rate')
    plt.legend(loc='best')
    pylab.savefig('output/sirf3d_vsAlpha'+timeNow())
    pylab.show()

    #fig=plt.figure()
    #CS = plt.contour(X, Y, Z)
    #plt.clabel(CS, inline=1, fontsize=10)
    #plt.title('attack rate')
    #plt.xlabel('alphas')
    #plt.ylabel('gammas')
    #CS = plt.contour(X, Y, Z)                                    
    ##CS = plt.contour(X, Y, Z, 6, colors='k', )
    #pylab.savefig('output/sirf3d_contours'+timeNow())
    #pylab.show()

    '''
    fig = plt.figure(figsize=plt.figaspect(2.))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    #ax = fig.add_subplot(2, 1, 1)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0, antialiased=False)
    plt.xlabel('alphas')
    plt.ylabel('gammas')
    pylab.savefig('output/sirf_3d'+timeNow())
    pylab.show()

    '''


def plot_synthetic():
    diamond = nx.DiGraph()
    #diamond.add_edges_from([(0,1),(0,2),(1,3),(2,3),(1,4),(2,4),(1,5),(2,5),(1,6),(2,6)])
    diamond.add_edges_from([(0,1),(0,2),(1,3),(2,3),(1,4),(2,4),(1,5),(2,5),])
    diamond.name = 'diamond '

    #configG = configuration_model()

    clique = nx.complete_graph(40)
    clique.name = 'K%d'%clique.number_of_nodes()
    
    #grid = nx.convert_node_labels_to_integers(nx.grid_2d_graph(10,5))
    #grid.name = 'grid10x5'

    #trirand = triangle_rand_g()
    graphs = [         diamond, ]
    #graphs = [         clique, ]
    names = ['SIR', '2FleeSIR']
    startNodes = [0]
    markers = ['r-^','b-o']
    linewidths = [1,2,3,4]
    pylab.figure()
    vsTauG = lambda G, tau, startNodes, transmissionModel: resMetrics.extentSIRSimLAdaptive(G, params={'tau':tau, 'vaccinated':[], 'transmissionModel':transmissionModel, 'weighted':False, 'minReplications':40, 'startNodes':startNodes, 'tolerance':0.005, 'K':dict.fromkeys(G.nodes(), 2)})[0]
    for j,G in enumerate(graphs):
        funcs = [lambda tau: vsTauG(G, tau, startNodes, 'SIR'), lambda tau: vsTauG(G, tau, startNodes, 'KFleeSIR')]
        taus = np.arange(0, 1.01, 0.05)
        pylab.hold(True)
        for i in xrange(2):
            #extents = np.array(map(funcs[i], taus))/G.number_of_nodes()
            extents = np.array(map(funcs[i], taus))
            pylab.plot(taus,extents,markers[i],linewidth=linewidths[i],label=G.name+'(%s)'%names[i])
            #pylab.plot(taus,map(lambda x:x*(clique.number_of_nodes()-1)+1, taus),label='line_fit')
    pylab.hold(False)
    pylab.legend(loc='best')
    pylab.xlabel(r'$\alpha$', fontsize=20)
    #pylab.ylabel('fraction of nodes infected', fontsize=20)
    pylab.ylabel('Number of nodes infected', fontsize=20)
    pylab.savefig('output/epidemics'+timeNow())
    pylab.show()

def predictors_of_flight():
#which factors predict flight response
    npr.seed(42)
    random.seed(42)
    import scipy.stats 

    #G = nx.complete_graph(40)
    #G.name = 'K%d'%G.number_of_nodes()
    nn  = 1000
    pER = 0.01
    G = nx.erdos_renyi_graph(nn, pER)
    G.name = 'ER(%d,%.3f)'%(nn, pER)
    #mAttach = 10
    #G = nx.barabasi_albert_graph(nn, mAttach)
    #G.name = 'BA(%d,%d)'%(nn, mAttach)
    
    #lam = 16.
    #lam = 10
    lam = 2.
    tau = 0.3

    Kvals = dict([(node,npr.poisson(lam)) for node in G])
    #Kvals = dict([(node,npr.binomial(G.number_of_nodes(), lam)) for node in G])
    
    extentDummy,statuses = resMetrics.extentKFleeSIRSim_helper(G=G, startN=0, vaccinated=[], weighted=False, verbose=True,
                                             params={'tau':tau, 'K':Kvals,})

    metrics = {'degree':nx.degree(G), 'betweenness':nx.betweenness_centrality(G)}
    metricData = []
    for node in G:
        nodeData = {'status':statuses[node][:3]}
        for metric in metrics:
            nodeData[metric] = metrics[metric][node]
        metricData.append(nodeData)

    for metric in metrics:
        print metric
        susceptible = [data[metric] for data in metricData if data['status']=='UIF']
        infected    = [data[metric] for data in metricData if data['status']=='INF']
        contrarian  = [data[metric] for data in metricData if data['status']=='CON']
        mn = min(susceptible + infected + contrarian)
        mx = max(susceptible + infected + contrarian)
        bins = np.linspace(mn, mx, num=10)
        #pdb.set_trace()
        pylab.figure()
        #pylab.title(metric)
        pylab.xlabel(metric, fontsize=20)
        pylab.hold(True)
        pylab.ylabel('Number of nodes infected', fontsize=20)
        '''
        if len(susceptible) > 0:
            #x,dummy=np.histogram(a=susceptible, bins=bins, normed=True)
            #pylab.bar(bins[:-1], x+1, color='g', label='susceptible', alpha=0.5)
            #(Pdb) x[0]*(dummy[0]-dummy[1])
            #-1.0
            dummy=pylab.hist(x=susceptible, bins=bins, color='g', label='susceptible', normed=0, alpha=0.5)
        if len(infected) > 0:
            dummy=pylab.hist(x=infected, bins=bins, color='r', label='infected', normed=0, alpha=0.5)
        if len(contrarian) > 0:
            dummy=pylab.hist(x=contrarian, bins=bins, color='b', label='contrarian', normed=0, alpha=0.5)
        '''
        pylab.hist(x=[susceptible, infected, contrarian], bins=bins, color=['g', 'r', 'b'], label=['susceptible', 'infected', 'contrarian'], normed=0, alpha=0.5, histtype='barstacked')

        KS = scipy.stats.ks_2samp(infected, contrarian)
        print 'KS: '
        print KS
        MW = scipy.stats.mannwhitneyu(infected, contrarian)
        print 'MW: '
        print MW

        pylab.legend(loc='best')
        pylab.hold(False)
        pylab.text(0.33, 1.075, '$G=%s, \\alpha=%.2f, \\lambda=%.1f, pval(inf,cont)MW=%.3f$'%(G.name,tau,lam,KS[1]), transform = pylab.gca().transAxes)
        pylab.savefig('output/covariates_G=%s_metric=%s_tau=%.2f_lam=%.1f__'%(G.name,metric,tau,lam)+timeNow()+'.pdf')

    npr.seed()
    random.seed()

def sampson_graph():
    import zipfile, cStringIO

    zf = zipfile.ZipFile('sampson_data.zip') # zipfile object
    e1=cStringIO.StringIO(zf.read('samplike1.txt')) # read info file
    e2=cStringIO.StringIO(zf.read('samplike2.txt')) # read info file
    e3=cStringIO.StringIO(zf.read('samplike3.txt')) # read info file
    G1=nx.read_edgelist(e1,delimiter='\t')
    G2=nx.read_edgelist(e2,delimiter='\t')
    G3=nx.read_edgelist(e3,delimiter='\t')

    return G1,G2,G3


def triangle_rand_g():
    #deg_tri=[[1,0],[1,0],[1,0],[2,0],[1,0],[2,1],[0,1],[0,1]]
    deg_tri=[[3,2],[1,2],[3,2],[2,1],[3,1],[2,1],[2,1],[2,2]]
    G = nx.random_clustered_graph(deg_tri)
    G=nx.Graph(G)
    G.remove_edges_from(G.selfloop_edges())

    G.name = 'random.tri'
    return G


def watch_sim():
    G = nx.complete_graph(40)
    G.name = 'K%d'%G.number_of_nodes()
    #G = nx.erdos_renyi_graph(20, 0.1)
    #G.name = 'ER%d'%G.number_of_nodes()
    
    #lam = 16.
    #lam = 10
    #lam = 1.
    lam = 3
    tau = 0.2

    Kvals = dict([(node,npr.poisson(lam)) for node in G])
    
    extentDummy,statuses = resMetrics.extentKFleeSIRSim_helper(G=G, startN=0, vaccinated=[], weighted=False, verbose=True,
                                             params={'tau':tau, 'K':Kvals,})

    pos = nx.drawing.layout.spring_layout(G, iterations=90)
    uninfected_color = 'green' #'black'
    infected_color   = 'red'
    contrarian_color = 'blue'

    uninfected_nodes = []
    infected_nodes = []
    contrarian_nodes = []
    for node in G:
        status = statuses[node]
        if 'UIF' in status:
            uninfected_nodes.append(node)
        elif 'INF' in status:
            infected_nodes.append(node)
        elif 'CON' in status:
            contrarian_nodes.append(node)
        else:
            raise ValueError, 'unknown status'


    reg_node_size = 200
    pylab.figure()
    pylab.hold(True)
    nx.draw_networkx_nodes(G, pos=pos, nodelist=uninfected_nodes, node_color=uninfected_color, node_size=reg_node_size, alpha=0.4)
    nx.draw_networkx_nodes(G, pos=pos, nodelist=infected_nodes,   node_color=infected_color,   node_size=reg_node_size, alpha=0.4)
    nx.draw_networkx_nodes(G, pos=pos, nodelist=contrarian_nodes, node_color=contrarian_color, node_size=reg_node_size, alpha=0.4)
    #nx.draw_networkx_edges(G, pos=pos, alpha=0.1)
    nx.draw_networkx_labels(G, pos=pos, labels=statuses)
    #pylab.grid(b=False)
    pylab.axis('off')
    #pylab.axes(frameon=False)
    pylab.text(-0.2, -0.2, 'nn=%d, tau=%.2f, lam=%.1f'%(G.number_of_nodes(),tau,lam))
    pylab.savefig('output/epidemic_final_tau=%.2f_lam=%.1f'%(tau,lam)+timeNow()+'.pdf')


if __name__ == '__main__':
    predictors_of_flight()
    #plot_sirf()

    #compare_stochastic_threshold()
    #plot_real()
    #plot_synthetic()
    #watch_sim()
