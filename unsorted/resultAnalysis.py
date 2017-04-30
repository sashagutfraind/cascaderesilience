'''
Analysis and plotting of the solutions

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


'''

import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
import matplotlib
#for running scripts silently:
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('PS')
if 'matplotlib.backends' not in sys.modules: matplotlib.use('pdf')
import matplotlib.pylab as pylab
#import pylab
import networkx as nx  #must be after pylab to stop loading pylab
import numpy as np
import numpy.random as npr
import random
import cPickle
import scipy
import time, os

import netDesign
import netDesignUnrooted
import netDesignCentrality
import netDesignGnp
import netDesignMulticell

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
#timeNow = lambda : str(time.localtime()[:-6]) + '_' + str(time.localtime()[-6:-3])

outDir  = 'output'
try:
    if not os.path.exists(outDir):
        os.mkdir(outDir)
except:
    print outDir +' cannot be created or already exists'

def analyze_centrality_scan(dossier=[]):

    tauVals      = []
    correlations = {}
    for paramName, paramVal, report in dossier['reports']:
        print '%s=%.3f'%(paramName,paramVal)
        print
        centDegList = []
        for soln in report['solutions']:
            optConfig = soln['optConfig']
            #Zvals.append(optConfig.cost)
            #Kvals.append(optConfig.configParams['k'])
            for net in optConfig.nets:
                centDegList += centrality_tuples(net)

        tauVals.append(paramVal)
        correlations[paramVal] = np.corrcoef(np.array(centDegList))[0,1]

    tauVals.sort()
    correlationsList = [correlations[tau] for tau in tauVals]
    #pylab.rc('text', usetex=True)
    pylab.plot(tauVals, correlationsList, 'r-o')
    #pylab.title(r'lambda=%2.1f'%(rat,))
    pylab.xlabel(r'$\tau$', fontsize=20)
    pylab.ylabel('Correlation', fontsize=20)
    #pylab.savefig('output/centDegCorr_vs_tau_'+timeNow()+'.eps')
    pylab.savefig('output/centDegCorr_vs_tau_'+timeNow()+'.pdf')


def analyze_centrality_scan_old(dossier=[]):
    correlations = []
    for paramName, paramVal, report in dossier['reports']:
        print '%s=%.3f'%(paramName,paramVal)
        print
        correlations = []
        for soln in report['solutions']:
            optConfig = soln['optConfig']
            #Zvals.append(optConfig.cost)
            #Kvals.append(optConfig.configParams['k'])
            for net in optConfig.nets:
                correlations += [centrality_correlation(net)]

        print 'Centrality-Degree Correlation:'
        print 'Average: %.3f'%np.average(correlations)
        print 'SD:      %.3f'%np.std(correlations)
        print '25-75%% range: %.3f-%.3f'%(correlations[len(correlations)/3], correlations[len(correlations)*2/3])

    #pylab.rc('text', usetex=True)
    #pylab.plot(results[:,0], results[:,1], '.')
    #pylab.title(r'lambda=%2.1f'%(rat,))
    #pylab.xlabel('Degree', fontsize=20)
    #pylab.ylabel('Betweenness', fontsize=20)

    

def analyze_plot_k_opt(dossier, plotK=True):
    #wishlist: this should be able to extract different features, not only k
    Xs    = dossier['Xs']
    Xname = dossier['Xname']
    Ys    = dossier['Ys']
    Yname = dossier['Yname']

    XsList = Xs.tolist()
    YsList = Ys.tolist() #needed for searching

    for paramName, paramVal, report in dossier['reports']:
        Xvals = []
        Yvals = []
        Zvals = []
        Kvals = []
        for soln in report['solutions']:
            optConfig = soln['optConfig']
            Xvals.append(optConfig.costParams[Xname])
            Yvals.append(optConfig.costParams[Yname])
            Zvals.append(optConfig.cost)
            Kvals.append(optConfig.configParams['k'])
        Zmesh = np.zeros(shape=(len(Ys),len(Xs)))
        if plotK:
            for i, x in enumerate(Xvals):
                y = Yvals[i]
                Zmesh[YsList.index(y),XsList.index(x)] = Kvals[i]
        else:
            for i, x in enumerate(Xvals):
                y = Yvals[i]
                Zmesh[YsList.index(y),XsList.index(x)] = Zvals[i]

        data={}
        data['Xs'] = Xs
        data['Ys'] = Ys
        #data['Zs'] = Zmesh
        #data['Xmesh'] = Xvals
        #data['Ymesh'] = Yvals
        data['Zmesh'] = Zmesh
        
        if plotK:
            #filename='output/contour__%s_%.2f_kvals_'%(paramName,paramVal)+timeNow()+'.eps'
            filename='output/contour__%s_%.2f_kvals_'%(paramName,paramVal)+timeNow()+'.pdf'
        else:
            #filename='output/contour__%s_%.2f_fitness_'%(paramName,paramVal)+timeNow()+'.eps'
            filename='output/contour__%s_%.2f_fitness_'%(paramName,paramVal)+timeNow()+'.pdf'

        #TODO: is it possible to get the contours to show logarithmically, or at least the k=2 contour
        contourPlot(data=data, caption='', Xname=Xname, Yname=Yname, filename=filename)
        #pylab.title('%s=%.2f'%(paramName,paramVal))

def analyze_tau_scan(fname = None, report = None, styles=None):
    #analyzes a single scan over tau values, producing a graphic
    import graphStats #here -  it might fail on the servers
    
    print 'analyzing...'
    reportData = openReport(fname=fname, report=report)
    report     = reportData['report']
    report_file= reportData['report_file']
    tauVals    = reportData['tauVals']
    design     = reportData['design']

    stats = {}
    stats['Xs']         = tauVals
    stats['Avg Degree'] = {}
    stats['fitness']    = {}
    dummyInstance       = design() 
    for param in dummyInstance.configParams:
        stats[param] = {}
    for tau in tauVals:
        stats['Avg Degree'][tau] = []
        stats['fitness'][tau]   = []    
        for param in dummyInstance.configParams:
            stats[param][tau] = []  

    paramSetNumbers = range(len(report['problemParamSets']))
    for paramSetNum in paramSetNumbers:
        optNetConfig = report['solutions'][paramSetNum]['optConfig']
        extraData    = report['solutions'][paramSetNum]['extraData']  #TODO: get the Z values
        tau          = report['problemParamSets'][paramSetNum]['costParams']['tau']
        nets         = optNetConfig.nets
        for net in nets:
            stats['Avg Degree'][tau] += [graphStats.avgDegree(net)]
            for param in dummyInstance.configParams:
                stats[param][tau] += [optNetConfig.configParams[param]]  #this doesn't give much variability...
        stats['fitness'][tau] = [netReport['cost']['fitness'] for netReport in optNetConfig.netReports]
        
    if fname!=None:
        data_fname = os.path.splitext(fname)[0]
    else:
        data_fname = 'output/tau_scan_report_'+timeNow()
    singleParamPlot(data=stats, paramName='Avg Degree', data_fname=data_fname) 
    singleParamPlot(data=stats, paramName='fitness', data_fname=data_fname) 
    for param in dummyInstance.configParams:
        singleParamPlot(data=stats, paramName=param, data_fname=data_fname, savefig=True)
    
    if design == netDesignUnrooted.unrootedDesign:
        multiParamPlot(Xs=tauVals, data=stats, seriesNames=['edgeFactorA', 'edgeFactorB', 'pairFactorA', 'pairFactorB'], data_fname=data_fname, styles=styles)
    elif design == netDesignCentrality.cdaDesign:
        multiParamPlot(Xs=tauVals, data=stats, seriesNames=['edgeFactorM', 'edgeFactorB'], data_fname=data_fname, styles=styles)
    
    if report_file != None:
        report_file.close()

def avgDegreePlot(data={}):
    singleParamPlot(data=data, paramName='avgDegree')


def centrality_correlation(G):
    results      = []
    centralities = nx.centrality.brandes_betweenness_centrality(G, normalized=False)
    for node in G:
        #if G.degree(node) == 1:  #if there is just one root node and the rest are leaves, get just one point
        #    continue
        #jitter for visualization:
        #results.append((G.degree(node)+npr.rand()*0.4, centralities[node]))        
        results.append((G.degree(node), centralities[node]))

    results = np.array(results)
    
    corr = scipy.corrcoef(results[:,0], results[:,1])[0,1]
    return corr

def centrality_tuples(G):
    results      = []
    centralities = nx.centrality.brandes_betweenness_centrality(G, normalized=False)
    for node in G:
        if G.degree(node) == 1:  
            continue
        #jitter for visualization:
        #results.append((G.degree(node)+npr.rand()*0.4, centralities[node]))        
        results.append((G.degree(node), centralities[node]))

    #results = np.array(results)
    #corr = scipy.corrcoef(results[:,0], results[:,1])[0,1]
    return results

def compare_reports(reportDatas, quantityNames='fitness', data_fname=None, plotmethod=pylab.plot, params=None):
    #produces a comparison of various series from various reports
    #   reportDatas   - a dictionary of multiple report data are given through a prior call to openReport
    #   quantityNames - either a single string, or a dictionary indexed by report key (the information to be plotted)

    #wishlist: improve the current architecture:
    #1. data should be collected on demand, rather than everytime

    import graphStats #here -  it might fail on the servers
    print 'comparing...'

    if data_fname == None:
        data_fname = 'comparison_' + timeNow()
    if type(quantityNames) is str:
        quantityNames = dict.fromkeys([name for name,data in reportDatas], quantityNames)

    stats = {}
    Xs    = {}
    for runName,reportData in reportDatas: 
        report     = reportData['report']
        report_file= reportData['report_file']
        tauVals    = reportData['tauVals']
        design     = reportData['design']

        reportStats = {}
        #TODO: perhaps just create the keys in quantity names? 
        reportStats['Xs']         = tauVals
        reportStats['fitness']    = {}
        reportStats['resilience'] = {}
        reportStats['efficiency'] = {}
        reportStats['Avg Degree'] = {}
        reportStats['Avg Degree+1'] = {}
        reportStats['log10 Avg Degree+1'] = {}
        #reportStats['Med Degree'] = {}
        dummyInstance       = design() 
        for param in dummyInstance.configParams:
            reportStats[param] = {}
        for tau in tauVals:
            reportStats['fitness'][tau]    = []    
            reportStats['resilience'][tau] = []    
            reportStats['efficiency'][tau] = []    
            for param in dummyInstance.configParams:
                reportStats[param][tau] = []  

        paramSetNumbers = range(len(report['problemParamSets']))
        for paramSetNum in paramSetNumbers:
            optNetConfig = report['solutions'][paramSetNum]['optConfig']
            tau          = report['problemParamSets'][paramSetNum]['costParams']['tau']
            nets         = optNetConfig.nets 
            for net in nets:
                for param in dummyInstance.configParams:
                    reportStats[param][tau] += [optNetConfig.configParams[param]]  #this doesn't give much variability...
            reportStats['fitness'][tau] = [netReport['cost']['fitness'] for netReport in optNetConfig.netReports]
            reportStats['resilience'][tau] = [netReport['cost']['resilience'] for netReport in optNetConfig.netReports]
            reportStats['efficiency'][tau] = [netReport['cost']['efficiency'] for netReport in optNetConfig.netReports]
            reportStats['Avg Degree'][tau] = [graphStats.avgDegree(net) for net in optNetConfig.nets]
            reportStats['Avg Degree+1'][tau] = [1+x for x in (reportStats['Avg Degree'][tau])]
            reportStats['log10 Avg Degree+1'][tau] = [np.log10(1+x) for x in (reportStats['Avg Degree'][tau])]
            #reportStats['Med Degree'][tau] = [graphStats.medDegree(net) for net in optNetConfig.nets]

        if report_file != None:
            report_file.close()

        #seriesName     = runName+ ': ' + quantityNames[runName] 
        seriesName     = runName 
        yseries        = reportStats[quantityNames[runName]]

        Xs[seriesName]    = tauVals
        stats[seriesName] = [np.average(yseries[tau]) for tau in tauVals]

    multiParamPlot(Xs=Xs, data=stats, seriesNames=[name for name,data in reportDatas], data_fname=data_fname, plotmethod=plotmethod, params=params)
    

def compare_reports_sensitivity(reportDatas, quantityNames, data_fname='comparison_' + timeNow(), plotmethod=pylab.plot, params=None):
    #produces a sensitivity analysis from various reports
    #   reportDatas   - a dictionary of multiple report data are given through a prior call to openReport
    #   quantityNames - either a single string, or a dictionary indexed by report key (thus we can compare either the same quantity or different quantities)

    #wishlist:
    #1. sometimes the parameter data is all we want - it should not compute all of the possible queries.  
    #2. allow the option of plotting the variance alongside the optimum 
    #3. performance hit: if we want sensitivity to several parameters, the network must be re-created for every call.

    print 'analyzing sensitivity...'

    if type(quantityNames) is str:
        quantityNames = dict.fromkeys([name for name,data in reportDatas], quantityNames)

    threshold = params.get('threshold', 0.95) #fraction of the optimum needed to count

    stats = {}
    Xs    = {}
    #overall pattern: pull all the information out, then store in stats just the data we want
    for runName,reportData in reportDatas: 
        report     = reportData['report']
        report_file= reportData['report_file']
        tauVals    = reportData['tauVals']
        design     = reportData['design']

        reportStats   = {}
        dummyInstance = design() 
        for quantity in dummyInstance.configParams:
            reportStats[quantity] = {}
        reportStats['resilience'] = {}
        reportStats['efficiency'] = {}
        reportStats['Avg Degree'] = {}
        reportStats['Avg Degree+1'] = {}
        reportStats['log10 Avg Degree+1'] = {}
        #reportStats['Med Degree'] = {}
        for tau in tauVals:
            for quantity in reportStats:
                reportStats[quantity][tau] = []  

        desiredDataName = quantityNames[runName]
        #this implementation is too dependent on the grid search algorithm
        #   * it should rely on report['solutions'][paramSetNum]['paramVals']
        for paramSetNum, paramSet in enumerate(report['problemParamSets']):
            optNetConfig     = report['solutions'][paramSetNum]['optConfig']
            fitnessThreshold = optNetConfig.cost * threshold
            tau              = report['problemParamSets'][paramSetNum]['costParams']['tau']

            paramScanValues = None
            if 'searchAlgParams' in paramSet and 'paramData' in paramSet['searchAlgParams']:
                paramScanAxes = paramSet['searchAlgParams']['paramData'] #list of points along each config parameter
            #elif 'paramVals' in report['solutions'][paramSetNum]:
            #    paramScanValues = report['solutions'][paramSetNum]['paramVals'] #dict, keyed same way as Zfull
            else:
                raise ValueError, 'Cannot find information about parameters'

            configParamList = compare_reports_sensitivity_helper(scores=report['solutions'][paramSetNum]['extraData']['Zfull'],
                                                                 paramScanValues=None,  #TODO: use this field, not the axes
                                                                 fitnessThreshold=fitnessThreshold,
                                                                 desiredDataName=desiredDataName,
                                                                 paramScanAxes=paramScanAxes,
                                                                 optNetConfig=optNetConfig)
            for config in configParamList:
                reportStats[desiredDataName][tau] += [config[desiredDataName]] 
            
        if report_file != None:
            report_file.close()


        #seriesName     = runName+ ': ' + quantityNames[runName] 
        seriesName     = runName 
        yseries        = reportStats[desiredDataName]

        Xs[seriesName]    = tauVals
        stats[seriesName] = []
        for tau in tauVals:
            yseries_tau = yseries[tau]
            if len(yseries[tau]) > 2:
                if len(yseries[tau]) < 10:
                    print 'STD for series %s, tau=%.3f contains just %d points within the threshold!'%(seriesName,tau,len(yseries[tau]))
                stats[seriesName] += [np.std(yseries[tau], ddof=1)]
            else:
                #raise ValueError, 'Cannot show STD for series %s, tau=%.3f because less than 2 points are in threshold!'%(seriesName,tau)
                print 'Inserting 0 for STD for series %s, tau=%.3f because only %d points are in threshold!'%(seriesName,tau,len(yseries[tau]))
                stats[seriesName] += [0.]
        
    multiParamPlot(Xs=Xs, data=stats, seriesNames=[name for name,data in reportDatas], data_fname=data_fname, plotmethod=plotmethod, params=params)


def compare_reports_sensitivity_helper(scores, paramScanValues, fitnessThreshold, desiredDataName, paramScanAxes=None, optNetConfig=None):
#1. generates list of configs (=a dict) where fitness passes the threshold
#2. creates a list of statistics, for each such config
    import graphStats
    configParamList = []

    #this identifies which parameter setting was made by the search algorithm
    #we could also just get all the parameters of the configuration
    if paramScanAxes!=None:
        paramScanValues = {}
        if len(paramScanAxes) == 1:
            XparamName = paramScanAxes[0][0]
            Xparams    = paramScanAxes[0][1]
            Xindices   = range(len(Xparams))
            for i in Xindices:
                paramScanValues[i] = {XparamName: Xparams[i]}
        elif len(paramScanAxes) == 2:
            XparamName = paramScanAxes[0][0]
            Xparams    = paramScanAxes[0][1]
            Xindices   = range(len(Xparams))
            YparamName = paramScanAxes[1][0]
            Yparams    = paramScanAxes[1][1]
            Yindices   = range(len(Yparams))
            for i in Xindices:
                for j in Yindices:
                    paramScanValues[(i,j)] = {XparamName: Xparams[i], YparamName: Yparams[j]}
        else:
            raise ValueError, 'No support for 3D grid search'

    availableParamNames = set([paramScanAxis[0] for paramScanAxis in paramScanAxes])
    availableParamNames.update(['resilience','efficiency'])
    needToReconstructNets = desiredDataName not in availableParamNames

    for paramIdx in paramScanValues:
        fitnesses = [net['cost']['fitness'] for net in scores[paramIdx]]
        if np.average(fitnesses) < fitnessThreshold:
            continue

        if needToReconstructNets:
            if optNetConfig != None: #recreate the net, and compute its avgDegree
                tmpInstance = optNetConfig.__class__(x=optNetConfig)  #__class__ will get the constructor 
                for param,val in paramScanValues[paramIdx].items():
                    tmpInstance.setParam(param,val)
                tmpInstance.buildNets(seeds=[netReport.get('seed',None) for netReport in scores[paramIdx]])
                avgDegree = np.average([graphStats.avgDegree(net) for net in tmpInstance.nets])
                paramScanValues[paramIdx]['Avg Degree'] = avgDegree
                paramScanValues[paramIdx]['Avg Degree+1'] = avgDegree+1.
                paramScanValues[paramIdx]['log10 Avg Degree+1'] = np.log10(1.+avgDegree)
                #medDegree = np.average([graphStats.medDegree(net) for net in tmpInstance.nets])
                #paramScanValues[paramIdx]['Med Degree'] = medDegree
            else:
                raise ValueError, 'To compute sensitivity, need to reconstruct net, but design and default param settings are not available'

        #this must be after the reconstruction of nets b/c it attempts to set parameters
        paramScanValues[paramIdx]['resilience'] = np.average([net['cost']['resilience'] for net in scores[paramIdx]])
        paramScanValues[paramIdx]['efficiency'] = np.average([net['cost']['efficiency'] for net in scores[paramIdx]])

        configParamList += [paramScanValues[paramIdx]]

    return configParamList


def contourPlot(data={}, caption='', Xname='', Yname='', filename=None):
    Xs = data['Xs'] #a single axis
    Ys = data['Ys'] #a single axis
    Zmesh = data['Zmesh']

    #Xmesh, Ymesh = np.meshgrid(Xs, Ys) #for some reason, this doesn't work if Xs and Ys have different lengths!

    import matplotlib.pyplot as plt
    plt.figure()
    CS = plt.contour(Xs, Ys, Zmesh)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title(caption)
    if Xname != 'tau':
        plt.xlabel(Xname, fontsize=20)
    else:
        plt.xlabel(r'$\tau$', fontsize=20)
    if Yname != 'tau':
        plt.ylabel(Yname, fontsize=20)
    else:
        plt.ylabel(r'$\tau$', fontsize=20)

    if filename==None:
        #filename = 'contour_'+timeNow()+'.eps'
        filename = 'contour_'+timeNow()+'.pdf'
    plt.savefig(filename)

    #see more features in:
    #http://matplotlib.sourceforge.net/plot_directive/mpl_examples/pylab_examples/contour_demo.py

def frontierPlot(scanDatas, fname = None, fullParamsSet = None, params=None):
    '''computes and prints the efficienct frontier from a set of simulations in report'''
    import graphStats #here -  it might fail on the servers
    import ensga
   
    class networkSample:
        def __init__(self, resilience, efficiency):
            self.dataVector = np.array([resilience, efficiency])
        def fitnessStats(self):
            return {'resilience':self.dataVector[0],'efficiency':self.dataVector[1]}
        def fitnessVector(self):
            return self.dataVector

    if params == None:
        params = {}
    print 'analyzing...'
    frontiers = []
    seriesNames = []
    for seriesName,scanData in scanDatas:
        samples = []
        for dataPt in scanData.values(): 
            resiliences  = [netReport['cost']['resilience'] for netReport in dataPt]
            efficiencies = [netReport['cost']['efficiency'] for netReport in dataPt]

            sample = networkSample(np.average(resiliences),np.average(efficiencies))
            samples.append(sample)

        arx = ensga.update(arx={}, samples=samples, 
                 saParams={'epsTolerance':params.get('epsTolerance', 0.03)})
        frontiers.append(arx)
        seriesNames.append(seriesName)
    
    if fname!=None:
        data_fname = os.path.splitext(fname)[0]
    else:
        data_fname = 'fronts_'+timeNow()

    write_frontiers(frontiers=frontiers, data_fname=data_fname, seriesNames=seriesNames, params=params)

def maxDegreePlot(data={}):
    singleParamPlot(data=data, paramName='maxDegree')

def multiParamPlot(Xs, data={}, xlab=r'$\tau$', seriesNames=[], data_fname='plot', plotmethod=pylab.plot, params=None):
    pylab.figure()
    pylab.hold(True)

    if not type(Xs) is dict:
        Xs_dict = {}
        for seriesName in data['jointSeriesNames']:
            Xs_dict[seriesName] = Xs
        Xs = Xs_dict

    if params==None:
        params = {}
    styles = params.get('styles', None)
    if styles==None:
        stylesAr  = ['k-|', 'b-s', 'g-o', 'm-^', 'c-D', 'r-v', 'y-p', 'b-<']
        #stylesAr = ['m-v', 'k--|', 'b--s', 'g-o', 'r-s', 'c--D']
        #stylesAr = ['k--|', 'b--s', 'g-o', 'm-s', 'c--D', 'r-v', 'y-s', 'b-o']
        #for longer ones:
        #stylesAr = ['k--', 'b--',  'g-',  'r-',  'c--', 'm--', 'y--', 'b-o']
        #sizesAr  = [9., 7., 5., 6., 3., 1., 8., 10.]
        styles = dict(zip(seriesNames, stylesAr[:len(seriesNames)]))
    sizes  = [1., 1., 1., 1., 1., 1., 1., 1.]

    pylab.xlabel(xlab, fontsize=30)
    ylab = params.get('title', None)
    logbase = params.get('logbase', 10)
    if ylab != None:
        pylab.ylabel(ylab, fontsize=20)

    minX = np.inf
    maxX = -np.inf
    minY = np.inf
    maxY = -np.inf
    for i, seriesName in enumerate(seriesNames):
        xseries = Xs[seriesName]
        yseries = data[seriesName]
        if plotmethod != pylab.semilogy:
            plotmethod(xseries, yseries, styles[seriesName], label=seriesName, linewidth=sizes[i])
        else:
            plotmethod(xseries, yseries, styles[seriesName], label=seriesName, linewidth=sizes[i], basey=logbase)

        minX = min(minX, min(xseries)-0.02*abs(max(xseries)-min(xseries)))
        maxX = max(maxX, 1.02*max(xseries))
        minY = min(minY, min(yseries)-0.02*abs(max(yseries)-min(yseries)))
        maxY = max(maxY, 1.02*max(yseries))
    minX = params.get('minX', minX)
    maxX = params.get('maxX', maxX)
    minY = params.get('minY', minY)
    maxY = params.get('maxY', maxY)
    
    pylab.axis([minX, maxX, minY, maxY])
    pylab.legend(loc='best')
    pylab.hold(False)
    #data_fname = data_fname+ '_multi' + '.pdf'
    data_fname = data_fname+ '_multi'
    pylab.savefig(data_fname)
    print 'Written: '+data_fname


def openReport(fname = None, report = None):
    #analyzes a single scan over tau values
    import graphStats #here -  it might fail on the servers

    if fname == None and report == None:
        raise ValueError, 'No data given for analysis'
    if report == None:
        report_file = open(fname, 'r')
        report      = cPickle.load(report_file)
    else:
        report_file = None

    tauVals = []
    for paramSet in report['problemParamSets']:
        tau          = paramSet['costParams']['tau']
        tauVals.append(tau)
    
    designName = report['problemParamSets'][0]['design'].split('.')
    design = getattr(globals()[designName[0]], designName[1])

    print 'Opened design report: ' + designName[1]

    reportData = {
                  'report':  report,
                  'report_file': report_file,
                  'tauVals': tauVals,
                  'design':  design
                 }

    return reportData


def singleParamPlot(data={}, paramName='', xlab=r'$\tau$', data_fname='', figure=None, savefig=True, updateAxes=True):
    #data[paramName] is keyed by paramVal, gives arrays to be plotted: 1->[1.1,0.9,1.0], 2->[3.1,2.9,3.11] etc
    if figure==None:
        figure = pylab.figure()
    pylab.hold(True)
    paramData = data[paramName]
    Xs = []
    Ys = []
    for xval in paramData:
        Xs = Xs + [xval] * len(paramData[xval]) 
        Ys = Ys + paramData[xval]
    pylab.plot(Xs, Ys, 'b.')

    Xs_avg = paramData.keys()
    Xs_avg.sort()
    Ys_avg = []
    for xval in Xs_avg:
        Ys_avg += [np.average(paramData[xval])]
    pylab.plot(Xs_avg, Ys_avg, 'b-', label=paramName)
    if updateAxes:
        pylab.axis([min(Xs)-0.02*abs(max(Xs)-min(Xs)), 1.02*max(Xs), min(Ys)-0.02*abs(max(Ys)-min(Ys)), 1.02*max(Ys)])

    pylab.xlabel(xlab, fontsize=20)
    pylab.legend()
    pylab.hold(False)
    if savefig:
        #pylab.savefig(data_fname + '__paramName_' + paramName + '.png')
        #pylab.savefig(data_fname + '__paramName_' + paramName + '.eps')
        pylab.savefig(data_fname + '__paramName_' + paramName + '.pdf')

def write_frontiers(frontiers, show=True, data_fname=None, figure=None, seriesNames=None, params=None):
    '''displays the list of efficient frontiers
    each front contains points with which are either coordinate pairs or has fitnessStats() method
    '''

    print 'plotting...'

    if type(params) is dict and 'styles' in params:
        styles = params['styles']
    else:
        styles = ['.', ',', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', '+', 'x', 'D', '|', '_']
        #colors  = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        #symbol  = markers[npr.randint(len(markers))] + colors[npr.randint(len(colors))]
    if seriesNames==None:
        seriesNames = [' ']*len(frontiers)

    if figure==None:
        figure = pylab.figure()
    pylab.hold(True)

    for i,front in enumerate(frontiers):
        if len(front) < 2:
            print 'Warning: Front %s has only %d points'%(seriesNames[i],len(front))
        pairs = []
        for resident in front:
            if hasattr(resident, 'fitnessStats'):
                fitness = resident.fitnessStats()
            else:
                fitness = resident
            pairs.append( (fitness['resilience'],fitness['efficiency']) )

        pairs.sort(key=lambda p:p[0])
        Xs = [pair[0] for pair in pairs]
        Ys = [pair[1] for pair in pairs]

        pylab.xlabel('Resilience', fontsize=20)
        pylab.ylabel('Efficiency', fontsize=20)
        if type(styles) is dict:
            pylab.plot(Xs, Ys, styles[seriesNames[i]], label=seriesNames[i])
        else:
            pylab.plot(Xs, Ys, styles[i], label=seriesNames[i])
        #time.sleep(2)
    pylab.legend(loc='best')
    pylab.title(params.get('title', ' '))
    if show:
        pylab.show()

    pylab.hold(False)
    if data_fname==None:
        data_fname = 'fronts_'+timeNow()
    pylab.savefig(data_fname+'.pdf')


if __name__ == '__main__':
    pass
