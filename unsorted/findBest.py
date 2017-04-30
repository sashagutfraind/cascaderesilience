'''
Search for the best network configuration

    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 


TODO:
    1. fail safe in case some simulations crash (failed file?)
    2. estimate for running time
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
import matplotlib.pylab as pylab
import networkx as nx  #must be after pylab, else gives problems with X
#import pylab
import time, os, sys
import cPickle
import pp
import taskDispatcher
import pdb

#import graphStats
import netDesign
import netDesignUnrooted
import netDesignCentrality
import netDesignGnp
import netDesignMulticell

import sa, search
import resultAnalysis


timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
#timeNow = lambda : str(time.localtime()[:-6]) + '_' + str(time.localtime()[-6:-3])



def compare_bifurcation_r_plot():
    try:
        dossier_file = open(r'output/scan_bifurcation_2009_07_14__16_44_31.pkl', 'r')
        dossier      = cPickle.load(dossier_file)
        resultAnalysis.analyze_plot_k_opt(dossier=dossier, plotK=True)
    except Exception, inst:
        print inst
        raise

def compare_centrality_plot():
    try:
        reports = []
        baseDir = r'sim-data\ns-tree-05-11-multiscan' 
        reports += [('tau', 0.00, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.000_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(9, 10, 24, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.05, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.050_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(9, 36, 49, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.10, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.100_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(10, 5, 28, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.15, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.150_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 13)_(2, 41, 28, 2, 133, 1).pkl')['report'])] 
        #reports += [('tau', 0.20, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.200_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(13, 39, 50, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.25, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.250_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(14, 2, 0, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.30, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.300_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(14, 28, 13, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.35, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.350_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(14, 58, 26, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.40, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.400_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(15, 37, 0, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.45, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.450_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(16, 26, 41, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.50, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.500_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(17, 35, 15, 0, 131, 1).pkl')['report'])] 
        #reports += [('tau', 0.55, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.550_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(19, 18, 31, 0, 131, 1).pkl')['report'])] 
        reports += [('tau', 0.60, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.600_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 11)_(21, 54, 24, 0, 131, 1).pkl')['report'])] 
        reports += [('tau', 0.65, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.650_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 13)_(6, 20, 42, 2, 133, 1).pkl')['report'])] 
        #reports += [('tau', 0.70, resultAnalysis.openReport(fname=baseDir+r'\')['report'])] 
        #reports += [('tau', 0.75, resultAnalysis.openReport(fname=baseDir+r'\')['report'])] 
        #reports += [('tau', 0.80, resultAnalysis.openReport(fname=baseDir+r'\')['report'])] 
        #reports += [('tau', 0.85, resultAnalysis.openReport(fname=baseDir+r'\')['report'])] 
        #reports += [('tau', 0.90, resultAnalysis.openReport(fname=baseDir+r'\multiScan_tau=0.900_design=netDesignCentrality.centralityTreeDesign_(2009, 5, 13)_(8, 21, 4, 2, 133, 1).pkl')['report'])] 
        #reports += [('tau', 0.95, resultAnalysis.openReport(fname=baseDir+r'\')['report'])] 
        #reports += [('tau', 1.00, resultAnalysis.openReport(fname=baseDir+r'\')['report'])] 

        dossier = {}
        dossier['reports'] = reports
        #dossier['Xname']   = Xname
        #dossier['Xs']      = Xs
        #dossier['Yname']   = Yname
        #dossier['Ys']      = Ys


        resultAnalysis.analyze_centrality_scan(dossier=dossier)
    except Exception, inst:
        print inst
        raise


def compare_g_vals():
    try:
        stylesAr  = ['k-|', 'b-s', 'g-o', 'm-^', 'c-D', 'r-v', 'y-p', 'b-<']
        for gval in [0.1, 10.]:
            if gval == 0.1:
                reportCavemen = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=0.1/tauScan_design=netDesignMulticell.cavemenDesign_(2009, 5, 13)_(22, 47, 56).pkl')
                reportConnStar= resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=0.1/tauScan_design=netDesignMulticell.connectedStarsDesign_(2009, 5, 13)_(7, 12, 47).pkl')
                reportClique  = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=0.1/tauScan_design=netDesignMulticell.cliqueDesign_(2009, 5, 13)_(7, 24, 9).pkl')
                reportCycle   = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=0.1/tauScan_design=netDesignMulticell.cycleDesign_(2009, 5, 13)_(7, 25, 6).pkl')      
                reportGnp     = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=0.1/tauScan_design=netDesignGnp.gnpDesign_(2009, 5, 13)_(7, 22, 49).pkl')
                reportStar    = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=0.1/tauScan_design=netDesignMulticell.starDesign_(2009, 5, 13)_(7, 26, 1).pkl')
                report911     = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=0.1/tauScan_design=netDesign.constantDesign_(2009, 5, 31)_(11, 1, 45)_911.pkl')
                reportFTP     = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=0.1/tauScan_design=netDesign.constantDesign_(2009, 5, 31)_(11, 3, 24)_FTP.pkl')
            elif gval == 10.:
                reportCavemen = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=10/')
                reportConnStar= resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=10/')
                reportClique  = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=10/')
                reportCycle   = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=10/')      
                reportGnp     = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=10/')
                reportStar    = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=10/')
                report911     = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=10/tauScan_design=netDesign.constantDesign_(2009, 5, 31)_(11, 6, 11)_911.pkl')
                reportFTP     = resultAnalysis.openReport(fname='sim-data/ns-scan-g/g=10/tauScan_design=netDesign.constantDesign_(2009, 5, 31)_(11, 7, 44)_FTP.pkl')
            else:
                return

            sources  = [
                          ('ConnCliques',   reportCavemen),
                          ('ConnStars', reportConnStar),
                          ('Cliques',   reportClique),
                          ('Cycles',    reportCycle),
                          ('ER',    reportGnp),
                          ('Stars',     reportStar),
                          ('9/11',      report911),
                          ('FTP',       reportFTP),
                          ]

            styles = dict(zip([name for name,rep in sources], stylesAr[:len(sources)]))  #used for all plots
            resultAnalysis.compare_reports(reportDatas=sources, quantityNames='fitness', data_fname='output/fitness_g-0_%d__'%int(1000*gval) + timeNow(), 
                                           params={'styles':styles, 'title':'Fitness'})
    
    except Exception, inst:
        print inst
        raise


def compare_gnp_g0_0_plot():
    reportGnp049    = resultAnalysis.openReport(fname='sim-data/ns-tau_gnp_g_eq_0-2009-08-19/tauScan_design=netDesignGnp.gnpDesign_2009_08_20__18_18_55.pkl')
    reportGnp051    = resultAnalysis.openReport(fname='sim-data/ns-tau_gnp_g_eq_0-2009-08-21/tauScan_design=netDesignGnp.gnpDesign_2009_08_21__15_37_48.pkl')
    sources    = [
                  ('ER_0.49',          reportGnp049),
                  ('ER_0.51',          reportGnp051),
                  ]


    stylesAr  = ['k-|', 'b-s', 'g-o', 'm-^', 'c-D', 'r-v', 'y-p', 'b-<']
    styles = dict(zip([name for name,rep in sources], stylesAr[:len(sources)]))  #used for all plots

    resultAnalysis.compare_reports(reportDatas=sources, quantityNames='fitness', data_fname='output/fitness_' + timeNow(), 
                                   params={'styles':styles, 'title':'Fitness'})

    resultAnalysis.compare_reports(reportDatas=sources, quantityNames='p', data_fname='output/p_multicell_' + timeNow(), 
                                                plotmethod=pylab.semilogy, 
                                                params={'styles':styles})

def compare_tau_scan_plot(r=0.51):
    try:
        if r==0.25:
            reportCavemen = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cavemenDesign_2010_02_24__15_57_45_transmutedR_0d25.pkl')
            reportConnStar= resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.connectedStarsDesign_2010_02_15__17_57_06_transmutedR_0d25.pkl')
            reportClique  = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cliqueDesign_2010_02_04__23_42_47_transmutedR_0d25.pkl')
            reportCycle   = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cycleDesign_2010_02_05__00_39_55_transmutedR_0d25.pkl')
            reportGnp     = resultAnalysis.openReport(fname='output/tauScan_design=netDesignGnp.gnpDesign_2010_02_05__00_32_43_transmutedR_0d25.pkl')
            reportStar    = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.starDesign_2010_02_05__01_07_53_transmutedR_0d25.pkl')
        elif r==0.49:
            reportCavemen = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cavemenDesign_2010_02_24__15_57_45_transmutedR_0d49.pkl')
            reportConnStar= resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.connectedStarsDesign_2010_02_15__17_57_06_transmutedR_0d49.pkl')
            reportClique  = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cliqueDesign_2010_02_04__23_42_47_transmutedR_0d49.pkl')
            reportCycle   = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cycleDesign_2010_02_05__00_39_55_transmutedR_0d49.pkl')
            reportGnp     = resultAnalysis.openReport(fname='output/tauScan_design=netDesignGnp.gnpDesign_2010_02_05__00_32_43_transmutedR_0d49.pkl')
            reportStar    = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.starDesign_2010_02_05__01_07_53_transmutedR_0d49.pkl')
        elif r==0.51:
            reportCavemen = resultAnalysis.openReport(fname='zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.cavemenDesign_2010_02_24__15_57_45.pkl')
            reportConnStar= resultAnalysis.openReport(fname='zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.connectedStarsDesign_2010_02_15__17_57_06.pkl')
            reportClique  = resultAnalysis.openReport(fname='zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.cliqueDesign_2010_02_04__23_42_47.pkl')
            reportCycle   = resultAnalysis.openReport(fname='zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.cycleDesign_2010_02_05__00_39_55.pkl')
            reportGnp     = resultAnalysis.openReport(fname='zoo/g1_r0d51/output/tauScan_design=netDesignGnp.gnpDesign_2010_02_05__00_32_43.pkl')
            reportStar    = resultAnalysis.openReport(fname='zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.starDesign_2010_02_05__01_07_53.pkl')
        elif r==0.75:
            reportCavemen = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cavemenDesign_2010_02_24__15_57_45_transmutedR_0d75.pkl')
            reportConnStar= resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.connectedStarsDesign_2010_02_15__17_57_06_transmutedR_0d75.pkl')
            reportClique  = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cliqueDesign_2010_02_04__23_42_47_transmutedR_0d75.pkl')
            reportCycle   = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cycleDesign_2010_02_05__00_39_55_transmutedR_0d75.pkl')
            reportGnp     = resultAnalysis.openReport(fname='output/tauScan_design=netDesignGnp.gnpDesign_2010_02_05__00_32_43_transmutedR_0d75.pkl')
            reportStar    = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.starDesign_2010_02_05__01_07_53_transmutedR_0d75.pkl')
        else:
            raise RuntimeError, 'Unknown value of r!'
       
        sources    = [
                      ('ConnCliques', reportCavemen),
                      ('ConnStars',   reportConnStar),
                      ('Cliques',     reportClique),
                      ('Cycles',      reportCycle),
                      ('ER',          reportGnp),
                      ('Stars',       reportStar),
                      ]
        stylesAr  = ['k-|', 'b-s', 'g-o', 'm-^', 'c-D', 'r-v', 'y-p', 'b-<']
        styles = dict(zip([name for name,rep in sources], stylesAr[:len(sources)]))  #used for all plots
        rStr   = 'R_0d%d_'%int(r*100) + timeNow()
        #resultAnalysis.compare_reports(reportDatas=sources, quantityNames='fitness', data_fname='output/fitness_'+rStr , 
        #                                           params={'styles':styles, 'title':'Fitness'})

        #resultAnalysis.compare_reports(reportDatas=sources, quantityNames='resilience', data_fname='output/resilience_' + rStr, 
        #                                   params={'styles':styles, 'title':'Resilience'})
        resultAnalysis.compare_reports_sensitivity(reportDatas=sources, quantityNames='resilience', data_fname='output/resilience_stdev_' + rStr, 
                                                    params={'threshold':.95, 'styles':styles, 'minY':0.0, 'maxY':1.0,})
        #resultAnalysis.compare_reports(reportDatas=sources, quantityNames='efficiency', data_fname='output/efficiency_' + rStr, 
        #                                   params={'styles':styles, 'title':'Efficiency'})
        resultAnalysis.compare_reports_sensitivity(reportDatas=sources, quantityNames='efficiency', data_fname='output/efficiency_stdev_' + rStr, 
                                                    params={'threshold':.95, 'styles':styles, 'minY':0.0, 'maxY':1.0,})

        #resultAnalysis.compare_reports(reportDatas=sources, quantityNames='Avg Degree', data_fname='output/avgDegree_'+rStr, 
        #                                   params={'styles':styles, 'title':'Avg Degree'})
        #resultAnalysis.compare_reports(reportDatas=sources, quantityNames='Avg Degree+1', data_fname='output/avgDegreePlus1_log_'+rStr, plotmethod=pylab.semilogy, 
        #                                   params={'styles':styles, 'title':'Avg Degree+1', 'minY': 0., 'maxY':185., })
        #resultAnalysis.compare_reports_sensitivity(reportDatas=sources, quantityNames='Avg Degree', data_fname='output/avgDegreePlus1_stdev_'+rStr, 
        #                                            params={'threshold':.95, 'title':'Avg Degree SD', styles':styles})
        #resultAnalysis.compare_reports_sensitivity(reportDatas=sources, quantityNames='Avg Degree+1', data_fname='output/avgDegreePlus1_stdev_log_'+rStr, plotmethod=pylab.semilogy,
        #                                    params={'threshold':.95, 'title':'Avg Degree+1 SD', 'styles':styles, 'minY': 0., 'maxY':185., })

        sources_k_cell    = [
                      ('ConnCliques', reportCavemen),
                      ('ConnStars',   reportConnStar),
                      ('Cliques',     reportClique),
                      ('Cycles',      reportCycle),
                      ('Stars',       reportStar),
                      ]
        #resultAnalysis.compare_reports(reportDatas=sources_k_cell, quantityNames='k', data_fname='output/k_multicell_' + rStr, plotmethod=pylab.semilogy, 
        #                                            params={'styles':styles, 'minY': 0., 'maxY': 181.})
        #resultAnalysis.compare_reports_sensitivity(reportDatas=sources_k_cell, quantityNames='k', data_fname='output/k_multicell_stdev_' + rStr, plotmethod=pylab.semilogy, 
        #                                            params={'threshold':.95, 'styles':styles, 'minY':0., 'maxY': 181.})

        sources_p    = [
                      ('ConnCliques', reportCavemen),
                      ('ConnStars',   reportConnStar),
                      ('ER',          reportGnp),
                      ]
        #resultAnalysis.compare_reports(reportDatas=sources_p, quantityNames='p', data_fname='output/p_' + rStr, 
        #                                            params={'styles':styles})
        resultAnalysis.compare_reports_sensitivity(reportDatas=sources_p, quantityNames='p', data_fname='output/p_stdev_' + rStr,
                                                    params={'threshold':.95, 'styles':styles, 'minY':0.0, 'maxY':1.0,})

       
    except Exception, inst:
        print inst
        raise


def compare_tau_scan_plot_empirical():
    r=0.51
    try:
        if r==0.49:
            #reportPower   = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_power.pkl')
            report11m     = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_11m_unwt.pkl')
            report911     = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_911.pkl')
            reportFTP     = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_ftp3c.pkl')
            reportEmail   = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_email.pkl') 
            reportNetsci  = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_netsci.pkl')
            reportAutoSys = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_autosys.pkl')
            reportGnutella= resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_gnutella.pkl')
            #reportConnStar= resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.connectedStarsDesign_2010_02_15__17_57_06_transmutedR_0d49.pkl')
        elif r==0.51:
            report11m     = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_11m_unwt.pkl')
            report911     = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_911.pkl')
            reportFTP     = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_ftp3c.pkl')
            reportEmail   = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_email.pkl') 
            reportNetsci  = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_netsci.pkl')
            reportAutoSys = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_autosys.pkl')
            reportGnutella= resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_gnutella.pkl')
            #reportConnStar= resultAnalysis.openReport(fname='zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.connectedStarsDesign_2010_02_15__17_57_06.pkl')
        else:
            raise RuntimeError, 'Unknown value of r!'
       
        sources    = []
        #if 'reportConnStar' in dir(): sources.append(('ConnStars',   reportConnStar))
        if 'report11m'      in dir(): sources.append(('11M',         report11m     ))
        if 'report911'      in dir(): sources.append(('9/11',        report911     ))
        if 'reportNetsci'   in dir(): sources.append(('CollabNet',   reportNetsci  ))
        if 'reportEmail'    in dir(): sources.append(('Email',       reportEmail   ))
        if 'reportFTP'      in dir(): sources.append(('FTP',         reportFTP     ))
        if 'reportGnutella' in dir(): sources.append(('Gnutella',    reportGnutella))
        if 'reportAutoSys'  in dir(): sources.append(('Internet AS', reportAutoSys ))
        #if 'reportPower'    in dir(): sources.append(('Power Grid',  reportPower   ))

        stylesAr  = ['b-s', 'k-|', 'g-o', 'm-^', 'c-D', 'r-v', 'y-p', 'b-<']
        styles = dict(zip([name for name,rep in sources], stylesAr[:len(sources)]))  #used for all plots
        rStr   = 'R_0d%d_'%int(r*100) + timeNow()

        resultAnalysis.compare_reports(reportDatas=sources, quantityNames='fitness', data_fname='output/fitness_emp_'+rStr, 
                                           params={'styles':styles, 'title':'Fitness'})#, 'minY':0.0, 'maxY':1.0})

        resultAnalysis.compare_reports(reportDatas=sources, quantityNames='resilience', data_fname='output/resilience_emp_'+rStr, 
                                           params={'styles':styles, 'title':'Resilience', 'minY':0.0, 'maxY':1.0})
        resultAnalysis.compare_reports(reportDatas=sources, quantityNames='efficiency', data_fname='output/efficiency_emp_'+rStr, 
                                           params={'styles':styles, 'title':'Efficiency', 'minY':0.0, 'maxY':1.0})

        
        #to get the correct BB use 'ps2eps -B -E filename.ps; epstopdf filename.eps'
    except Exception, inst:
        print inst
        raise
       
def compare_tau_scan_plot_empirical_binary():
    r=0.51
    try:
        if r==0.49:
            report11m_unwt = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_11m_unwt.pkl')
            report11m_wted = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_11m_wted.pkl')
            report911_unwt  = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_911.pkl')
            report911_wted  = resultAnalysis.openReport(fname='sim-data/g1_r0d49_empirical/tauScan_911wted.pkl')
        elif r==0.51:
            report11m_unwt = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_11m_unwt.pkl')
            report11m_wted = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_11m_wted.pkl')
            report911_unwt  = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_911.pkl')
            report911_wted  = resultAnalysis.openReport(fname='sim-data/g1_r0d51_empirical/tauScan_911wted.pkl')
        else:
            raise RuntimeError, 'Unknown value of r!'
       
        sources    = []
        if 'report11m_unwt' in dir(): sources.append(('11M  (binary)',   report11m_unwt     ))
        if 'report11m_wted' in dir(): sources.append(('11M  (weighted)', report11m_wted     ))
        if 'report911_unwt' in dir(): sources.append(('9/11 (binary)',   report911_unwt     ))
        if 'report911_wted' in dir(): sources.append(('9/11 (weighted)', report911_wted     ))

        stylesAr  = ['bs', 'b-s', 'ro', 'r-o', 'cD', 'c-D', 'yp', 'y-p']
        styles = dict(zip([name for name,rep in sources], stylesAr[:len(sources)]))  #used for all plots
        rStr   = 'R_0d%d_'%int(r*100) + timeNow()

        resultAnalysis.compare_reports(reportDatas=sources, quantityNames='fitness', data_fname='output/fitness_emp_'+rStr, 
                                           params={'styles':styles, 'title':'Fitness', 'minY':0.0, 'maxY':1.0})

        resultAnalysis.compare_reports(reportDatas=sources, quantityNames='resilience', data_fname='output/resilience_emp_'+rStr, 
                                           params={'styles':styles, 'title':'Resilience', 'minY':0.0, 'maxY':1.0})
        resultAnalysis.compare_reports(reportDatas=sources, quantityNames='efficiency', data_fname='output/efficiency_emp_'+rStr, 
                                           params={'styles':styles, 'title':'Efficiency', 'minY':0.0, 'maxY':1.0})

       
    except Exception, inst:
        print inst
        raise


def computeNets(problemParamSets = [], pathDir = 'output', simName = '', job_server=None, batchSize=40, debugThreads=False):
#finds the best network using a search algorithm for a given parameter set
#threaded version
    report = {}
    report['problemParamSets'] = problemParamSets
    report['timeNow']          = timeNow()
    
    def helperFn(args):
        paramSet = args['paramSet']
        return computeNetInstance(paramSet)


    storage_lock = taskDispatcher.thread.allocate_lock()
    threads = {}
    results = {}
    tasks = range(len(problemParamSets))
    while len(tasks) > 0:
        problemParamSetsBatch = tasks[:batchSize]
        for paramIdx in problemParamSetsBatch:
           th = taskDispatcher.DictionaryStorageThread(fn=helperFn, args={'paramSet':problemParamSets[paramIdx]}, storage=results, 
                                                       result_key=paramIdx, storage_lock=storage_lock, debug=debugThreads)
           threads[paramIdx] = th
        #separate loop to ensure they all have the lock before start
        for th_key in threads:
           threads[th_key].start()

        while threads != {}:
            paramIdx, th = threads.items()[0]
            if th.isAlive():
               time.sleep(1)
            else:
               threads.pop(paramIdx)

        tasks = tasks[batchSize:]
    #print storage_lock
    del storage_lock

    #prune to enable pickling
    for paramSet in problemParamSets:
        paramSet['searchAlg'] = str(paramSet['searchAlg'])
        paramSet['design']    = str(paramSet['design'])
        paramSet['costParams']['costFn'] = str(paramSet['costParams']['costFn']) 
        costParams = paramSet['costParams']
        if 'job_server' in costParams: #cannot be pickled
            costParams.pop('job_server')
    #copy results and prune
    report['solutions'] = [None]*len(problemParamSets)
    for paramIdx, paramSet in enumerate(problemParamSets):
        paramSetResults = results[paramIdx]
        report['solutions'][paramIdx] = paramSetResults

        optConfig = paramSetResults['optConfig']
        optConfig.costParams.pop('job_server')
        if 'job_server' in costParams: #cannot be pickled
            report[paramIdx].costParams.pop('job_server')
        
    #print report.keys()
    #print report

    filename   = simName + '_' + timeNow() + '.pkl'
    safePickle(dirname=pathDir, filename=filename, data=report)

    return report

def computeNetsDispatch(problemParamSets = [], pathDir = 'output', simName = '', job_server=None):
    """
    runs each of the problems using an indirect calling approach:
    1. each paramSet is written to an input .pkl file.
    2. using PP or PBS, taskDispatcher.py is invoked at each of the worker nodes and given the input
       (to prevent blocking of this controller we invoke using threads)
       the input .pkl contains the name of the output file we expect to receive
    3. taskDispatcher.py calls the appropriate algorithm
    4. taskDispatcher.py stores the output in a .pkl, which we collect in an infinite loop

    """
    import os, cPickle, subprocess
    from threading import Thread

    inputFilenames = []
    outputFilenames = []
    resultsIndices = {}
    print 'Preparing input parameter files...'
    for paramSetNum, paramSet in enumerate(problemParamSets):
        algParam = paramSet.copy()
        algParam['costParams']['job_server'] = None #otherwise cannot be passed through pp
        
        timeRec = timeNow()
        inputFilename  = pathDir + os.sep + 'input_'+timeRec+'_paramSetNum'+str(paramSetNum)+'.pkl'
        outputFilename = pathDir + os.sep + 'output_'+timeRec+'_paramSetNum'+str(paramSetNum)+'.pkl'
        inputFilenames.append(inputFilename) 
        outputFilenames.append(outputFilename) 
        resultsIndices[outputFilename] = paramSetNum

        executionParams = {}
        executionParams['fnFullName']     = 'findBest.computeNetInstance' 
        executionParams['paramSet']       = algParam
        executionParams['outputFilename'] = outputFilename

        inputFile = open(inputFilename, 'w')  
        cPickle.dump(executionParams, inputFile)
        inputFile.close()

    sys.stdout.flush()
    def jobCaller(params):
        try:
            outNerr = open(params[1], 'w')
            retCode = subprocess.call([params[0]], stdout=outNerr, stderr=outNerr, shell=True)
        finally:
            outNerr.close()
        return retCode

    print 'Dispatching calls...'
    jobs = []
    for paramSetNum, paramSet in enumerate(problemParamSets):
        cmdLineParams = (['python -s taskDispatcher.py --Pkl '+inputFilenames[paramSetNum], outputFilenames[paramSetNum]+'.txt'],)
        job = job_server.submit(func=jobCaller, args=cmdLineParams, depfuncs=(), modules=('subprocess',))
        jobs.append(job)
        #for PDB-system:
        #PID=os.spawnl(mode=P_NOWAIT, file='python', args='tastDispatcher.py --arg1=val1 --input=inputParams.pkl >& out.txt')

    #threads are useful only for pp, because job() must be called, yet it would block us from retrieving the outputs
    #  they chaperone the jobs
    class jobThread(Thread):
       def __init__ (self,job,jobNum):
         Thread.__init__(self)
         self.job=job
         self.jobNum=jobNum
       def run(self):
         print 'Starting job '+str(self.jobNum)
         retCode = self.job()
         if retCode != 0:
            print 'Computation did not return 0 for job #'+str(self.jobNum)
            print retCode
   
    sys.stdout.flush()
    jobThreads = []
    try:
        for jobNum, job in enumerate(jobs):
            th = jobThread(job=job,jobNum=jobNum)
            th.start()
            jobThreads.append(th)
    except:
        print 'Error in starting threads!'
        #th.join(timeout=1)

    for inputFilename in inputFilenames:
        if not os.path.exists(inputFilename):
            raise RuntimeError, 'Input file was not created! '%inputFilename

    results = [None]*len(problemParamSets)
    
    print
    print 'Waiting for calls to complete...'
    print
    soughtOutputs = set(outputFilenames)
    soughtTxtputs = set([name+'.txt' for name in outputFilenames])
    reportProgress = True
    failedOutputs = []
    failedReturn = lambda outData: outData['optConfig'].cost == None or not np.isfinite(outData['optConfig'].cost)
    startTime = time.time()
    while len(soughtOutputs) > 0:
        time.sleep(1)
        for txtputFilename in list(soughtTxtputs):
            if os.path.exists(txtputFilename):
                soughtTxtputs.remove(txtputFilename)
        if reportProgress and len(soughtTxtputs) == 0:
            print 
            print 'Progress: All scripts started!'
            print
            reportProgress = False

        for outputFilename in list(soughtOutputs):
            if os.path.exists(outputFilename+'.Failed'):
                print 'Error: output not available due to error '+outputFilename
            if not os.path.exists(outputFilename+'.Done'):
                continue
            with open(outputFilename, 'r') as f:
                outData = cPickle.load(f)
                if failedReturn(outData):
                    fFailed = open(outputFilename+'.Failed', 'w')
                    fFailed.close()
                    failedOutputs.append(outputFilename)
            results[resultsIndices[outputFilename]] = outData
            soughtOutputs.remove(outputFilename)
            #os.remove(inputFilenames[resultsIndices[outputFilename]])
            #os.remove(outputFilename)
            print 'Successfully pulled ' + outputFilename
        elapsed = int(time.time()-startTime)
        if elapsed%60==0:
            print 'Time: %dmin. Num remaining: %d'%(elapsed/60,len(soughtOutputs)) 
            sys.stdout.flush()

    print 'The following output jobs failed:'
    for repName in failedOutputs:
        print repName
        print 
    if len(failedOutputs)==0: print 'NONE'

    for th in jobThreads:
        th.join(timeout=0)

    sys.stdout.flush()

    report = {}
    report['problemParamSets'] = problemParamSets
    report['timeNow']          = timeNow()

    report['solutions'] = [None]*len(problemParamSets)
    for paramIdx, paramSet in enumerate(problemParamSets):
        #to enable pickling
        paramSet['searchAlg'] = str(paramSet['searchAlg'])
        paramSet['design']    = str(paramSet['design'])
        paramSet['costParams']['costFn'] = str(paramSet['costParams']['costFn']) 

        paramSetResults = results[paramIdx]
        optConfig = paramSetResults['optConfig']
        optConfig.costParams['job_server'] = str(optConfig.costParams.get('job_server', None))

        report['solutions'][paramIdx] = paramSetResults


    filename   = simName + '_' + timeNow() + '.pkl'
    safePickle(dirname=pathDir, filename=filename, data=report)

    return report

def computeNetInstance(paramSet):
    #optimization of network for a specified parameter set
    #may be called directed or from taskDispatcher.callMethod()
    tau     = paramSet['costParams'].get('tau', 0.5)
    print 'The tau=' + str(tau) + ' case:'

    alg     = paramSet.get('searchAlg', sa.runSA)

    ###
    optNetConfig = alg(paramSet)
    ###

    if 'debug' in paramSet['searchAlgParams'].get('outputLevel', set()):
        print 'The optimal configuration is:'
        if type(optNetConfig) is not dict:
            print optNetConfig
        else:
            print optNetConfig['optConfig']  #just the opt, without extraData
        print 'Done.'
    else:
        print 'Done.'
        print 

    return optNetConfig


def frontier_plot(targetTaus=None):
    reportCavemen = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cavemenDesign_2010_02_24__15_57_45_transmutedR_0d49.pkl')
    reportConnStar= resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.connectedStarsDesign_2010_02_15__17_57_06_transmutedR_0d49.pkl')
    reportClique  = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cliqueDesign_2010_02_04__23_42_47_transmutedR_0d49.pkl')
    reportCycle   = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.cycleDesign_2010_02_05__00_39_55_transmutedR_0d49.pkl')
    reportGnp     = resultAnalysis.openReport(fname='output/tauScan_design=netDesignGnp.gnpDesign_2010_02_05__00_32_43_transmutedR_0d49.pkl')
    reportStar    = resultAnalysis.openReport(fname='output/tauScan_design=netDesignMulticell.starDesign_2010_02_05__01_07_53_transmutedR_0d49.pkl')

    if targetTaus == None:
        targetTaus = [0.9]
    sources    = [
                  ('ConnCliques', reportCavemen['report']),
                  ('ConnStars',   reportConnStar['report']),
                  ('Cliques',     reportClique['report']),
                  ('Cycles',      reportCycle['report']),
                  ('ER',          reportGnp['report']),
                  ('Stars',       reportStar['report']),
                  ]

    #sources = [
    #            ('Cycles', reportCycle['report']),  #finds that all resiliences close to 1
    #            ('ER', reportGnp['report']),
    #            ('Stars', reportStar['report'])
    #          ]

    for targetTau in targetTaus:
        scanDatas = []
        for reportName,report in sources:
            paramSetNumbers = range(len(report['problemParamSets']))
            for paramSetNum in paramSetNumbers:
                #searchAlgParams = report['searchAlgParams'][paramSetNum]['paramData']
                problemParams   = report['problemParamSets'][paramSetNum]  #check: does this include the variable params
                scanData    = report['solutions'][paramSetNum]['extraData']['Zfull']
                tau         = problemParams['costParams']['tau']
                if tau == targetTau:
                    scanDatas.append((reportName,scanData))
                    break
            
        stylesAr  = ['k-|', 'b-s', 'g-o', 'm-^', 'c-D', 'r-v', 'y-p', 'b-<']
        styles = dict(zip([name for name,rep in sources], stylesAr[:len(sources)]))  #used for all plots
        
        resultAnalysis.frontierPlot(scanDatas=scanDatas, fname='output/frontiers_tau0d%d_'%(int(targetTau*100),)+timeNow(), fullParamsSet=problemParams, params={'styles':styles, 
                                'epsTolerance':0.03, 
                                #'title':r'$\tau$=%.3f'%targetTau
                                })
        

def produceNetFigure(netDesignClass, s, label):
    npr.seed(s)
    net = netDesignClass.buildNet()
    pylab.figure()
    nx.draw(net)
    pylab.savefig('output/'+label)
    return net


def run_bifurcation_r_stars():
    job_server = taskDispatcher.create_job_server(simulatedPP=True) #simulated
    #job_server = taskDispatcher.create_job_server(simulatedPP=False)
    
    try:
        dossierStars = scan_bifurcation_r(job_server=job_server, design=netDesignMulticell.starDesign)
        resultAnalysis.analyze_plot_k_opt(dossier=dossierStars)
    except Exception, inst:
        print inst
        raise

    if job_server != None:
        job_server.destroy()


def run_centrality():
    #job_server = taskDispatcher.create_job_server(simulatedPP=True) #simulated
    job_server = taskDispatcher.create_job_server(simulatedPP=False)
    
    try:
        dossierCent = scan_centrality(job_server=job_server, design=netDesignCentrality.centralityTreeDesign)
        resultAnalysis.analyze_centrality_scan(dossier=dossierCent)
    except Exception, inst:
        print inst
        raise

    if job_server != None:
        job_server.destroy()

def run_empirical():
    job_server = taskDispatcher.create_job_server(simulatedPP=True) #simulated
    #job_server = taskDispatcher.create_job_server(simulatedPP=False)
    r=0.51  #warning: not used by all sims below
    '''
    try:
        G = nx.Graph(nx.read_gml('data/G911_unwted.gml'))
        G = nx.convert_node_labels_to_integers(G, first_label=0)
        reportG = scan_tau(job_server=job_server, design=netDesign.constantDesign, fixedParams={'G':G}, costParams={'weighted':False, 'r':r})
        resultAnalysis.analyze_tau_scan(report=reportG, fname='output/tauScan_911unwted'+timeNow())
    except Exception, inst:
        print inst
        raise
    try:
        G = nx.Graph(nx.read_gml('data/G911_wted.gml'))
        G = nx.convert_node_labels_to_integers(G, first_label=0)
        reportG = scan_tau(job_server=job_server, design=netDesign.constantDesign, fixedParams={'G':G}, costParams={'weighted':True, 'r':r})
        resultAnalysis.analyze_tau_scan(report=reportG, fname='output/tauScan_911wted'+timeNow())
    except Exception, inst:
        print inst
        raise
    '''
    try:
        G = nx.Graph(nx.read_gml('data/G11m_unwted.gml'))
        G = nx.convert_node_labels_to_integers(G, first_label=0)
        reportG = scan_tau(job_server=job_server, design=netDesign.constantDesign, fixedParams={'G':G}, costParams={'weighted':False, 'r':r})
        resultAnalysis.analyze_tau_scan(report=reportG, fname='output/tauScan_11Munwted'+timeNow())
    except Exception, inst:
        print inst
        raise
    try:
        G = nx.Graph(nx.read_gml('data/G11m_wted.gml'))
        G = nx.convert_node_labels_to_integers(G, first_label=0)
        reportG = scan_tau(job_server=job_server, design=netDesign.constantDesign, fixedParams={'G':G}, costParams={'weighted':True, 'r':r})
        resultAnalysis.analyze_tau_scan(report=reportG, fname='output/tauScan_11Mwted'+timeNow())
    except Exception, inst:
        print inst
        raise

    '''
    try:
        Gftp3c = nx.read_edgelist('data/ftp3c.graph', edgetype=int)
        Gftp3c = nx.convert_node_labels_to_integers(Gftp3c, first_label=0)
        reportGftp3c = scan_tau(job_server=job_server, design=netDesign.constantDesign, fixedParams={'G':Gftp3c})
        #resultAnalysis.analyze_tau_scan(report=reportGftp3c, fname='output/tauScan_ftp3c_'+timeNow())
    except Exception, inst:
        print inst
        raise

    try:
        Gftp = nx.read_edgelist('data/ftp.graph', edgetype=int)
        Gftp = nx.convert_node_labels_to_integers(Gftp, first_label=0)
        reportGftp = scan_tau(job_server=job_server, design=netDesign.constantDesign, GfixedParams={'G':Gftp})
        resultAnalysis.analyze_tau_scan(report=reportGftp, fname='output/tauScan_ftp_'+timeNow())
    except Exception, inst:
        print inst
        raise

    try:
        Gemail = nx.read_edgelist('data/arenas_email.edges', nodetype=int)
        Gemail = nx.convert_node_labels_to_integers(Gemail, first_label=0)
        reportGemail = scan_tau(job_server=job_server, design=netDesign.constantDesign, GfixedParams={'G':Gemail})
        resultAnalysis.analyze_tau_scan(report=reportGemail, fname='output/tauScan_email'+timeNow())
    except Exception, inst:
        print inst
        raise
    
    try:
        Gnetsci = nx.Graph(nx.read_gml('data/newman06_netscience.gml'))
        Gnetsci = nx.convert_node_labels_to_integers(Gnetsci, first_label=0)
        reportGnetsci = scan_tau(job_server=job_server, design=netDesign.constantDesign, GfixedParams={'G':Gnetsci})
        resultAnalysis.analyze_tau_scan(report=reportGnetsci, fname='output/tauScan_netScience'+timeNow())
    except Exception, inst:
        print inst
        raise
    
    try:
        Gpower = nx.Graph(nx.read_gml('data/watts_strogatz98_power.gml'))
        reportGpower = scan_tau(job_server=job_server, design=netDesign.constantDesign, GfixedParams={'G':Gpower})
        resultAnalysis.analyze_tau_scan(report=reportGpower, fname='output/tauScan_powerGrid'+timeNow())
    except Exception, inst:
        print inst
        raise

    try:
        GautoSys = nx.read_edgelist('data/as-caida20071105.txt', edgetype=int).to_undirected()
        GautoSys = nx.convert_node_labels_to_integers(GautoSys, first_label=0)
        reportAutoSys = scan_tau(job_server=job_server, design=netDesign.constantDesign, GfixedParams={'G':GautoSys})
        #resultAnalysis.analyze_tau_scan(report=reportAutoSys, fname='output/tauScan_autosys'+timeNow())
    except Exception, inst:
        print inst
        raise
    
    try:
        Ggnutella = nx.read_edgelist('data/p2p-Gnutella08.txt', edgetype=int).to_undirected()
        Ggnutella = nx.convert_node_labels_to_integers(Ggnutella, first_label=0)
        reportGnutella = scan_tau(job_server=job_server, design=netDesign.constantDesign, GfixedParams={'G':Ggnutella})
        #resultAnalysis.analyze_tau_scan(report=rreportGnutella, fname='output/tauScan_gnutella'+timeNow())
    except Exception, inst:
        print inst
        raise
    ''' 

    if job_server != None:
        job_server.destroy()


def run_tau_scan():
    simulatedPP = False

    job_server = taskDispatcher.create_job_server(simulatedPP=simulatedPP)
    try:
        reportClique = scan_tau(job_server=job_server, design=netDesignMulticell.cliqueDesign)
        resultAnalysis.analyze_tau_scan(report=reportClique)
        #resultAnalysis.analyze_tau_scan(fname='output/tauScan_design=netDesignMulticell.cliqueDesign_(2009, 4, 4)_(1, 33, 8).pkl')
    except Exception, inst:
        print inst
        raise
    job_server.destroy()

    job_server = taskDispatcher.create_job_server(simulatedPP=simulatedPP)
    try:
        reportGnp = scan_tau(job_server=job_server, design=netDesignGnp.gnpDesign)
        resultAnalysis.analyze_tau_scan(report=reportGnp)
        #resultAnalysis.analyze_tau_scan(fname='output/tauScan_design=netDesignGnp.gnpDesign_(2009, 4, 4)_(1, 34, 34).pkl')
    except Exception, inst:
        print inst
        raise
    job_server.destroy()

    job_server = taskDispatcher.create_job_server(simulatedPP=simulatedPP)
    try:
        reportCycle = scan_tau(job_server=job_server, design=netDesignMulticell.cycleDesign)
        #resultAnalysis.analyze_tau_scan(report=reportCycle)
        #resultAnalysis.analyze_tau_scan(fname='output/tauScan_design=netDesignMulticell.cycleDesign_(2009, 4, 4)_(1, 39, 30).pkl')
    except Exception, inst:
        print inst
        raise
    job_server.destroy()
    
    job_server = taskDispatcher.create_job_server(simulatedPP=simulatedPP)
    try:
        reportStar = scan_tau(job_server=job_server, design=netDesignMulticell.starDesign)
        resultAnalysis.analyze_tau_scan(report=reportStar)
        #resultAnalysis.analyze_tau_scan(fname='output/tauScan_design=netDesignMulticell.starDesign_(2009, 4, 4)_(1, 37, 30).pkl')
    except Exception, inst:
        print inst
        raise
    job_server.destroy()
   
    job_server = taskDispatcher.create_job_server(simulatedPP=simulatedPP)
    try:
        reportConStar = scan_tau(job_server=job_server, design=netDesignMulticell.connectedStarsDesign)
        #resultAnalysis.analyze_tau_scan(report=reportConStar)
    except Exception, inst:
        print inst
        raise
    job_server.destroy()

    job_server = taskDispatcher.create_job_server(simulatedPP=simulatedPP)
    try:
        reportCavemen = scan_tau(job_server=job_server, design=netDesignMulticell.cavemenDesign)
        #resultAnalysis.analyze_tau_scan(report=reportCavemen)
        #resultAnalysis.analyze_tau_scan(fname='data/netstorm-mu=0.499/tauScan_design=netDesignMulticell.cavemenDesign_(2009, 4, 6)_(16, 48, 4).pkl')
    except Exception, inst:
        print inst
        raise
    job_server.destroy()


def run_fitness_analytic():
    job_server = taskDispatcher.create_job_server(simulatedPP=True) #simulated
    #job_server = taskDispatcher.create_job_server(simulatedPP=False) #real

    try:
        dossier = scan_analytic(job_server=job_server, design=netDesignMulticell.starDesign)
        resultAnalysis.analyze_plot_k_opt(dossier=dossier, plotK=False)
    except Exception, inst:
        print inst
        raise


def safePickle(dirname, filename, data):
    try:
        if not os.path.exists(dirname):
            os.mkdir(dirname)
    except:
        print dirname +' cannot be created or already exists'

    try:
        fullpath   = dirname + os.sep + filename
        outputFile = open(fullpath, 'wb')
        cPickle.dump(data, outputFile)
        outputFile.close()
        print 
        print 'Pickle: ' + fullpath + ' written!'
    except Exception, inst:
        print 'Unable to pickle...'
        print inst

def scan_2d(job_server=None, design=netDesignMulticell.starDesign, scanFixedName='g', scanFixedVal=1.0, Xname='tau', Xs=None, Yname='r', Ys=None, analytic=True, configParams=None, costParams=None, fixedParams=None):
    print 
    print '%s: %.3f'%(scanFixedName,scanFixedVal)
    print 
    print 'Scan domain:'
    print 'X: %s'%Xname
    print Xs
    print 'Y: %s'%Yname
    print Ys

    if configParams==None:
        configParams = {}
    if fixedParams==None:
        fixedParams = {}
    if costParams==None:
        costParams = {}

    problemParamSets = []
    defaultParamSet  = {}
    defaultParamSet['design']           = design
    
    #WARNING: these are set for ALL networks
    defaultParamSet['fixedParams']      = {
                                           'nn':          180,
                                           'outputLevel': set(['debug']),
                                           }
    defaultParamSet['costParams']       = {'analytic':   analytic,
                                           'costFn':     '',
                                           'job_server':  job_server,
                                           #'curWorkDir':  os.getcwd(),
                                           }
    defaultParamSet['costParams'][scanFixedName] = scanFixedVal

    defaultParamSet['searchAlg']        = search.gridSearch
    defaultParamSet['searchAlgParams']  = {
                                           'simulatedSearch': False,
                                           'outputLevel': set(['debug', 'normal', 'extraData']),
                                           }
    if not analytic:
        defaultParamSet['fixedParams']['numNets']        = 10
        defaultParamSet['costParams']['costFn']          = netDesign.defaultSIRcost
        defaultParamSet['costParams']['job_server']      = job_server
        defaultParamSet['costParams']['minReplications'] = 40
        defaultParamSet['costParams']['tolerance']       = 0.01

    if design == netDesignMulticell.cliqueDesign \
         or design == netDesignMulticell.cycleDesign \
         or design == netDesignMulticell.starDesign:
        defaultParamSet['searchAlgParams']['paramData'] = [('k',np.array(range(0,181,1)))]
    elif design == netDesignCentrality.centralityTreeDesign:
        defaultParamSet['searchAlgParams']['paramData'] = [('q',np.arange(0, 3.01, 0.10)), ('b',np.array(np.arange(1.0, 20.01, 0.5).tolist() + [20, 33, 50, 100]))]
        #defaultParamSet['searchAlgParams']['paramData'] = [('q',np.arange(0, 2.01, 0.3)), ('b',np.array([10, 100]))]
        defaultParamSet['fixedParams']['numNets']        = 1 #all are the same
    else:
        raise ValueError, 'No analytic expression known for design: %s'%str(design)

    for k in configParams:
        defaultParamSet['configParams'][k] = configParams[k]

    for k in costParams:
        defaultParamSet['costParams'][k] = costParams[k]

    for k in fixedParams:
        defaultParamSet['fixedParams'][k] = fixedParams[k]

    #defaultParamSet['searchAlg']        = sa.runSA 
    print 'defaultParamSet:'
    print defaultParamSet
    print

    for x in Xs:
        for y in Ys:
            paramSet = {}
            for key in defaultParamSet:
                if hasattr(defaultParamSet[key], 'copy'):
                    paramSet[key] = defaultParamSet[key].copy()
                else:
                    paramSet[key] = defaultParamSet[key]
            paramSet['costParams'][Xname] = x
            paramSet['costParams'][Yname] = y
            problemParamSets.append(paramSet)

    debugThreads = False
    if debugThreads:
        print 
        print 'Warning: debugging threads!'
        print
    else:
        print 'Running multiple threads!'        
    simName = 'multiScan__%s_%.3f_design_'%(scanFixedName,scanFixedVal)+str(design)
    simName = simName.replace('.', '_') #the dot is trouble for latex
    #report  = computeNetsDispatch( 
    #                     problemParamSets = problemParamSets, 
    #                     pathDir          = 'output', 
    #                     simName          = simName, 
    #                     job_server       = job_server)

    report  = computeNets( 
                         problemParamSets = problemParamSets, 
                         pathDir          = 'output', 
                        simName          = simName, 
                         job_server       = job_server, 
                         batchSize        = 40, 
                         debugThreads     = debugThreads
                        )
    return report


def scan_analytic(job_server=None, design=netDesignMulticell.starDesign):
    print 
    print 'Initializing scan over g values:'
    print 'Design: ' + str(design)
    print 

    scanFixedName = 'g'
    #gVals = [0.0, 0.1, 0.5, 1.0, 2.0, 10.0]
    #gVals = [0.1, 10.0]
    gVals = [0.5, 2.0]

    Xname = 'tau'
    Xs    = np.arange(0, 1.01, 0.05)

    Yname = 'r'
    Ys    = np.arange(0, 1.01, 0.05)  
    
    reports = []
    for g in gVals:
        report = scan_2d(job_server=job_server, design=design, Xname=Xname, Xs=Xs, Yname=Yname, Ys=Ys, scanFixedName=scanFixedName, scanFixedVal=g, analytic=True)

        reports += [('g', g, report)]

    dossier = {}
    dossier['reports'] = reports
    dossier['Xname']   = Xname
    dossier['Xs']      = Xs
    dossier['Yname']   = Yname
    dossier['Ys']      = Ys

    return dossier

def scan_bifurcation_r(job_server=None, design=netDesignMulticell.starDesign):
    print 
    print 'Initializing scan...'
    print 'Design: ' + str(design)
    print 

    scanFixedName = 'g'
    gVals = np.array([1.0])

    Xname = 'tau'  #inside resultAnalysis it gets converted to \tau
    Xs    = np.arange(0.0, 1.01, 0.05) 
    
    Yname = 'r'
    Ys    = np.arange(0, 1.01, 0.05)
    
    reports = []
    for gVal in gVals:
        report = scan_2d(job_server=job_server, design=design, Xname=Xname, Xs=Xs, Yname=Yname, Ys=Ys, scanFixedName=scanFixedName, scanFixedVal=gVal, analytic=True)

        reports += [(scanFixedName, gVal, report)]

    dossier = {}
    dossier['reports'] = reports
    dossier['Xname']   = Xname
    dossier['Xs']      = Xs
    dossier['Yname']   = Yname
    dossier['Ys']      = Ys

    filename = 'scan_bifurcation_' + timeNow() + '.pkl'
    
    safePickle(dirname='output', filename=filename, data=dossier)

    return dossier

def scan_centrality(job_server=None, design=netDesignCentrality.centralityTreeDesign):
    print 
    print 'Initializing scan...'
    print 'Design: ' + str(design)
    print 

    scanFixedName = 'tau'
    tauVals    = np.arange(0.0, 1.01, 0.05) 

    Xname = 'g'
    Xs    = [0.0, 0.1, 0.5, 1.0, 2.0, 10.0]  

    Yname = 'r'
    Ys    = np.arange(0, 1.01, 0.05)
    
    reports = []
    for tau in tauVals:
        report = scan_2d(job_server=job_server, design=design, Xname=Xname, Xs=Xs, Yname=Yname, Ys=Ys, scanFixedName=scanFixedName, scanFixedVal=tau, analytic=False)

        reports += [(scanFixedName, tau, report)]

    dossier = {}
    dossier['reports'] = reports
    #dossier['Xname']   = Xname
    #dossier['Xs']      = Xs
    #dossier['Yname']   = Yname
    #dossier['Ys']      = Ys

    filename = 'scan_centrality_' + timeNow() + '.pkl'
    
    safePickle(dirname='output', filename=filename, data=dossier)

    return dossier


def scan_tau(job_server=None, design=netDesignGnp.gnpDesign, configParams=None, costParams=None, fixedParams=None):
    print 
    print 'Initializing scan over tau values:'
    print 'Design: ' + str(design)
    print 
    tauVals = np.arange(0, 1.01, 0.05)  
    
    print 'taus:'
    print tauVals

    if configParams==None:
        configParams = {}
    if fixedParams==None:
        fixedParams = {}
    if costParams==None:
        costParams = {}

    problemParamSets = []
    defaultParamSet  = {}
    defaultParamSet['design']           = design
    
    if design != netDesign.constantDesign:
        defaultParamSet['fixedParams']      = {
                                               'nn':          180, 
                                               'numNets':     10, 
                                               }
    else:
        G = fixedParams['G']
        defaultParamSet['fixedParams']= {}
        print 'Running on a fixed graph:'
        print '%d nodes and %d edges'%(G.number_of_nodes(),G.number_of_edges())

    defaultParamSet['fixedParams']['outputLevel'] = set(['debug', 'normal', 'extraData'])

    defaultParamSet['costParams']       = {'costFn':          netDesign.defaultSIRcost,
                                           'job_server':      job_server,
                                           #'curWorkDir':     os.getcwd(),
                                           'minReplications': 40,
                                           'tolerance':       0.005,
                                           'weighted':       False, #is the network weighted

                                           #'g':           0.0,
                                           'g':           1.0, 
                                           #'g':           10.0,

                                           #'r':          0.25} 
                                           'r':          0.49}  
                                           #'r':          0.51}
                                           #'r':          0.75}
    #defaultParamSet['searchAlg']        = search.nelderMead
    defaultParamSet['searchAlg']        = search.gridSearch
    defaultParamSet['searchAlgParams']  = {
                                           'simulatedSearch': False,
                                           #'runs': 2,
                                           'outputLevel': set(['debug', 'extraData']),
                                           }
    if design == netDesignMulticell.cavemenDesign \
            or design == netDesignMulticell.connectedStarsDesign:
        defaultParamSet['searchAlgParams']['paramData'] = [('p',np.array(np.arange(0.002, 0.034, 0.002).tolist() + np.arange(0, 1.01, 0.05).tolist()) ), 
                                                           ('k',np.arange(0,181,1))]
    elif design == netDesignMulticell.cliqueDesign \
         or design == netDesignMulticell.cycleDesign \
         or design == netDesignMulticell.starDesign: #k in star is already corrected for.  k=1 means just the center
        defaultParamSet['searchAlgParams']['paramData'] = [('k',np.arange(0,181,1))]
        defaultParamSet['fixedParams']['numNets']       = 1
    elif design == netDesignGnp.gnpDesign:
        defaultParamSet['searchAlgParams']['paramData'] = [('p',np.array(np.arange(0.002, 0.034, 0.002).tolist() + np.arange(0, .51, 0.05).tolist() + [1.0]) ), 
                                                          ]
    elif design == netDesignCentrality.centralityTreeDesign:
        defaultParamSet['searchAlgParams']['paramData'] = [('q',np.arange(0, 3.01, 0.10)), ('b',np.array(np.arange(1.0, 20.01, 0.5).tolist() + [50]))]
    elif design == netDesign.constantDesign:
        defaultParamSet['fixedParams']['G']       = G
        defaultParamSet['fixedParams']['nn']      = G.number_of_nodes()
        defaultParamSet['fixedParams']['numNets'] = 1

        defaultParamSet['searchAlgParams']['paramData'] = [('dummy', np.array([0]))]
    #elif design == netDesignCentrality.cdaDesign: #WARNING: fixedParams could lead to side effects with class's defaults
    #    defaultParamSet['searchAlgParams']['paramData'] = [('edgeFactorM',np.arange(-15.0, 15.01, 0.5)), ('edgeFactorB',np.arange(-15.0, 15.01, 0.5))]
    #    defaultParamSet['fixedParams']['numUpdateCycles']       = 10.
    #    defaultParamSet['fixedParams']['transitivityBias']      = 0.5
    #    defaultParamSet['fixedParams']['updatesPerCentrality']  = 5. 


    for k in configParams:
        defaultParamSet['configParams'][k] = configParams[k]

    for k in costParams:
        defaultParamSet['costParams'][k] = costParams[k]

    for k in fixedParams:
        defaultParamSet['fixedParams'][k] = fixedParams[k]

    #defaultParamSet['searchAlg']        = sa.runSA 
    print 'defaultParamSet:'
    print defaultParamSet

    for tau in tauVals:
        paramSet = {}
        for key in defaultParamSet:
            if hasattr(defaultParamSet[key], 'copy'):
                paramSet[key] = defaultParamSet[key].copy()
            else:
                paramSet[key] = defaultParamSet[key]
        paramSet['costParams']['tau'] = tau
        problemParamSets.append(paramSet)

    debugThreads = False
    #debugThreads = True
    simName      = 'tauScan_design='+str(design)
    report = computeNetsDispatch( 
                         problemParamSets = problemParamSets, 
                         pathDir          = 'output', 
                         simName          = simName, 
                         job_server       = job_server, 
                         )
    #report = computeNets( 
    #                     problemParamSets = problemParamSets, 
    #                     pathDir          = 'output', 
    #                     simName          = simName, 
    #                     job_server       = job_server, 
    #                     batchSize        = 40, 
    #                     debugThreads     = debugThreads
    #                    )
    return report

def transmute_r(targetR=0.49):
    fnames = ('zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.cycleDesign_2010_02_05__00_39_55.pkl',
              'zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.cavemenDesign_2010_02_24__15_57_45.pkl', 
              'zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.connectedStarsDesign_2010_02_15__17_57_06.pkl',
              'zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.cliqueDesign_2010_02_04__23_42_47.pkl',
              'zoo/g1_r0d51/output/tauScan_design=netDesignGnp.gnpDesign_2010_02_05__00_32_43.pkl',
              'zoo/g1_r0d51/output/tauScan_design=netDesignMulticell.starDesign_2010_02_05__01_07_53.pkl',
             )
    reports = [fname for fname in fnames]
    
    fitness = lambda data, targetR: (1-targetR)*data['efficiency'] + targetR*data['resilience'] 
    for fname in reports:
        with open(fname, 'r') as report_file:
            report = cPickle.load(report_file)
            paramSetNumbers = range(len(report['problemParamSets']))
            for paramSetNum in paramSetNumbers:
                problemParams   = report['problemParamSets'][paramSetNum]
                problemParams['costParams']['r'] = targetR

                solutions = report['solutions'][paramSetNum]
                optConfig = solutions['optConfig']
                extraData = solutions['extraData']
                Zfull     = extraData['Zfull']
                Zs        = extraData['Zs']
                paramVals = extraData['paramVals']
                bestConfigNum = -1
                bestConfigFitness = -np.inf
                for configNum in paramVals:
                    netFitnesses = []
                    for netNum,network in enumerate(Zfull[configNum]):
                        fit = fitness(data=Zfull[configNum][netNum]['cost'], targetR=targetR)
                        Zfull[configNum][netNum]['cost']['fitness'] = fit
                        netFitnesses.append(fit)
                    Zs[configNum] = np.average(netFitnesses)
                    if Zs[configNum] > bestConfigFitness:
                        bestConfigFitness = Zs[configNum]
                        bestConfigNum     = configNum
                for paramName,val in paramVals[bestConfigNum].items():
                    optConfig.setParam(paramName,val)
                optConfig.setParam('r', targetR)
                if Zfull[bestConfigNum][0].get('seed', False):
                    seeds = [report['seed'] for netReport in Zfull[bestConfigNum]]
                else:
                    seeds = None
                optConfig.costParams['job_server']=None
                optConfig.evaluate(seeds=seeds)
                print

            safePickle(dirname='output', filename=os.path.splitext(os.path.split(fname)[1])[0]+'_transmutedR_0d%d'%int(targetR*100)+'.pkl', data=report)
            del report
    return
        

if __name__ == '__main__':
    #transmute_r(targetR=0.25)
    #transmute_r(targetR=0.49)
    #transmute_r(targetR=0.75)

    #run_centrality()
    run_empirical()
    #run_tau_scan()
    #run_bifurcation_r_stars()
    #run_fitness_analytic()

    #frontier_plot(targetTaus=[0.3,0.4,0.9])

    #compare_gnp_g0_0_plot()
    #compare_tau_scan_plot(r=0.25)
    #compare_tau_scan_plot(r=0.49)
    #compare_tau_scan_plot(r=0.51)
    #compare_tau_scan_plot(r=0.75)
    #compare_tau_scan_plot_empirical()
    #compare_tau_scan_plot_empirical_binary()
    #compare_g_vals()
    #compare_centrality_plot()
    #compare_bifurcation_r_plot()

