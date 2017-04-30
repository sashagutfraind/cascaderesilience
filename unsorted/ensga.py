'''
eps-Non-dominated Sorting Genetic Algorithm for Multi-Objective Optimization

    ref: Laumanns, Thiele, Deb and Zitzler

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
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('pdf')
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('PS')
#import matplotlib.pylab as pylab
#import matplotlib.pyplot as plt
import networkx as nx
import pylab  #FIXME
import time, os, sys
import cPickle

import netDesign
import netDesignMulticell
import netDesignGnp
import resultAnalysis


def ensgaRun(inputParams):
    algParams                = inputParams.get('searchAlgParams', {})
    saParams = {}
    saParams['design']       = algParams['design']
    saParams['fixedParams']  = inputParams.get('fixedParams', {})
    saParams['costParams']   = inputParams.get('costParams', {})
    saParams['numSamples']   = algParams.get('numSamples', 20)
    saParams['epsTolerance'] = algParams.get('epsTolerance', 0.05)
    saParams['mutateRate']   = algParams.get('mutateRate', 0.2)
    saParams['outputLevel']  = algParams.get('outputLevel', set())
    saParams['maxIter']      = algParams.get('maxIter', 10E6) 
    saParams['callback']     = algParams.get('callback', lambda x: None) 

    outputLevel     = saParams['outputLevel']
    maxIter         = saParams['maxIter'] 
    callback        = saParams['callback']
    #simulatedSearch = algParams.get('simulatedSearch',   False)
    #paramData       = algParams.get('paramData', None)  #[(name1, range1), (name2, range2), ..]
    #if paramData == None:
    #    raise ValueError, 'Empty search problem!'

    if 'normal' in outputLevel:
        print 'Starting eNSGA Search Algorithm:'

    try:
        arx = {}
        extraData = {}
        parents = None
        iter = 0
        while (not terminate(arx)) and iter < maxIter:
            children = generate(arx=arx, parents=parents, saParams=saParams)
            arx      = update(arx=arx, samples=children, saParams=saParams)
            parents  = children
            iter += 1
            callback(parents, arx)

    except RuntimeError, e:
        print e
        print 'Error: cannot complete run!'
        print 'Likely b/c network couldn\'t be completed.'

    if 'normal' in outputLevel:
        print 'Optimal set: ', arx
    sys.stdout.flush()
    
    if 'extraData' in outputLevel:
        return {'archive':arx, 'extraData':extraData}
    else:
        return arx


def generate(arx, parents, saParams):
    design     = saParams['design']
    numSamples = saParams['numSamples']
    mutateRate = saParams['mutateRate']
    outputLevel= saParams['outputLevel']
    
    children = []
    if parents == None:
        for i in xrange(numSamples):
            child   = design(costParams=saParams['costParams'], fixedParams=saParams['fixedParams'])
            for m in xrange(40):
                child.mutate()
            try:
                child.evaluate()
            except Exception, e:
                if 'debug' in outputLevel:
                    print 'construction/evaluation failed:'
                    print e
                raise
            children.append(child)
        return children

    numParents = len(parents)
    if numParents == 0: raise RuntimeError, 'Empty parents array!'
    numTrials  = 10*numSamples
    for i in xrange(numTrials):
        if len(children) >= numSamples: 
            break

        p1num = npr.randint(numParents)
        p2num = npr.randint(numParents)
        parent1 = parents[p1num]
        parent2 = parents[p2num]

        child = parent1.cross(parent2)
        child.mutate(rate=mutateRate)
        try:
            child.evaluate()
        except Exception, e:
            if 'debug' in outputLevel:
                print 'construction/evaluation failed:'
                print e
            raise
        children.append(child)
    return children

 
def sampleBox(sample, epsTolerance):
    v = sample.fitnessVector() #assumes returns vector >= 0
    v = np.log(v)/np.log(1+epsTolerance) #if epsTolerance=0 returns -inf, no crash
    return np.floor(v)
   
def terminate(arx):
#a stopping condition based on the properties of the existing solutions
    return False

def update(arx, samples, saParams):
    '''updates the archive with new samples
       wishlist: check archive consistency, e.g. not more than two samples in the same box, no dominated pts
    '''
    epsTolerance = saParams['epsTolerance']

    for sample in samples:
        sBox = sampleBox(sample, epsTolerance)

        cDominatedResidents = set()
        fDominatedResident  = None #should be only one
        sampleInNewBox      = True
        sampleIsDominated   = False

        #option: maybe rather than using arx[resident]->resBox, we use arx[resBox]->resident
        for resident in arx:
            resBox = arx[resident]
            diffBox = sBox - resBox  #+ means sample is better than resident
            if np.abs(diffBox).sum() == 0.:
                sampleInNewBox = False
                if (sample.fitnessVector()-resident.fitnessVector()).min() > 0:
                    fDominatedResident = resident
                    break
            elif diffBox.min() >= 0:
                cDominatedResidents.add(resident)
            elif diffBox.max() <= 0:
                sampleIsDominated = True
                break
        if sampleIsDominated:
            if len(cDominatedResidents) != 0 or fDominatedResident != None:
                raise RuntimeError, 'Something went wrong: arx inconsistent?'
            continue
        if len(cDominatedResidents) != 0:
            for resident in cDominatedResidents:
                arx.pop(resident)
            arx[sample] = sBox
        if fDominatedResident != None:
            arx.pop(fDominatedResident)
            arx[sample] = sBox
        if sampleInNewBox:
            arx[sample] = sBox
    return arx            



def _test():
    inputParams = {
                   'costParams': {'analytic':True, 'g':1.0, 'tau':0.5},
                   'searchAlgParams': {'mutateRate':0.5,
                                       'design': netDesignMulticell.starDesign,
                                       'maxIter': 100,
                                       }
                   }

    arx = ensgaRun(inputParams)
    resultAnalysis.write_frontiers([arx])

def _test2():
    writeFrontWrap = lambda parents, arx: writeFront(arx=arx, show=False)
    inputParams = {
                   'costParams': {#'analytic':True, 
                                 'costFn':netDesign.defaultSIRcost,
                                 'minReplications':40,
                                 'tolerance':0.05,
                                 'g':1.0, 
                                 'tau':0.1,
                                 'r':0.5 #a stub
                                 },
                   'fixedParams': {'nn':180, 'numNets':5, 'outputLevel':set()},
                   'searchAlgParams': {'mutateRate':0.5,
                                       #'design': netDesignMulticell.starDesign,
                                       #'design': netDesignGnp.gnpDesign,
                                       'design': netDesignMulticell.connectedStarsDesign,
                                       'maxIter': 30,
                                       'epsTolerance':0.15,
                                       'callback':writeFrontWrap,
                                       }
                   }

    #pylab.show()
    pylab.figure()
    pylab.show() #THIS must be killed
    pylab.figure()
    pylab.hold(True)
    arx = ensgaRun(inputParams)

    resultAnalysis.write_frontiers([arx])

if __name__ == '__main__':
    _test2()

