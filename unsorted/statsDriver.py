"""
OBSOLETE

Run analyze data from simulations

TODO : implement!


    Copyright (C) 2007-2010 by  
    Sasha Gutfraind <ag362@cornell.edu> 
    Distributed under the terms of the GNU Lesser General Public License 
    http://www.gnu.org/copyleft/lesser.html 



Usage help
----------

python statsDriver.py -h

Release Notes
-------------

""" 
#    Copyright (C) 2007 by  
#    Sasha Gutfraind <ag362@cornell.edu> 
#    Distributed under the terms of the GNU Lesser General Public License 
#    http://www.gnu.org/copyleft/lesser.html 
# 
import sys
sys.path.append('/home/gfriend/lib/python')
sys.path.append('/homes/gfriend/lib/python')
import numpy
import matplotlib
#for running scripts silently:
#matplotlib.use('PS')
import matplotlib.pylab as pylab
#import pylab
import scipy
import re, os, sys

try:
    import rpy
except ImportError:
    print 'Couldn\'t import RPy. Statistical analysis functionality is deactivated.'

def usage():
    print "Script for [running the ABM and] analyzing the data"
    print "Allowed options are:"
    #print "[-b  <Batch file>]  Batch file to execute. Default: analysis1.pf"
    print "[-d  <Dir>]         Perform analysis on Dir containing the simulation output files. "
    print "[-a  <y|n|o>]       [y=perform batch and then analyze (default),]"
    print "                    [n=just run and not analyze,] o=analyze Only."
    print "-h                  Displays this"
    print "eg."
    print r"To run and analyze a batch file"
    print r"python abmdriver.py -b analysis1.pf"
    print
    print r"To analyze without running the batch file(data generated in the past)"
    print r"python abmdriver.py -b analysis1.pf -a o"


def initialize():
    import getopt
    opts, args = getopt.getopt(sys.argv[1:], "b:d:a:h", [""])

    batchFile = ''
    dataDir   = ''
    analysisMode = 'runANDanalyze'
    
    for o, a in opts:
       if o in ("-h", "--help"):
          usage()
          sys.exit(0)
       if o in ("-b"):
          batchFile = a
       if o in ("-d"):
          dataDir = a
          analysisMode = 'analyzeOnly'
       if o in ("-a"):
          if (a == 'y') or (a == 'Y'):
              analysisMode = 'runANDanalyze'
          elif (a == 'n') or (a == 'N'):
              analysisMode = 'runANDstop'
          elif (a == 'o') or (a == 'O'):
              analysisMode = 'analyzeOnly'
          else:
              raise ValueError, 'Couldn\'t parse analysis mode'

    if dataDir == '':
        raise IOError, 'No data dir specified'
        sys.exit(1)
    
    if not os.path.exists(dataDir):
        print 'outputDirectory "%s" doesn\'t exist' %dataDir
        sys.exit(1)
    
    return (batchFile, dataDir, analysisMode) 


def runRegression(statName, bestValues, seriesNames, f):
    f.write('\nRegression for ' + statName + '\n')

    xvals = numpy.array([], dtype=numpy.double)
    yvals = numpy.array([], dtype=numpy.double)
    f.write(statName + '\t' + 'Mean\n')
    for i, name in enumerate(seriesNames):
        Ys = numpy.array(bestValues[i], dtype=numpy.double)
        try:
            Xs = float(seriesNames[i]) * numpy.ones(Ys.size, dtype=numpy.double)
        except:
            Xs = float(i)              * numpy.ones(Ys.size, dtype=numpy.double)
        xvals = numpy.concatenate((xvals, Xs))
        yvals = numpy.concatenate((yvals, Ys))

        f.write(seriesNames[i] + '\t' + str(Ys.mean()) + '\n')

    data={'x':xvals, 'y':yvals}
    lm_res = rpy.r.lm(rpy.r('y~x'), data)
    lm_sum = rpy.r.summary_lm(lm_res)
    lm_aov = rpy.r.summary_aov(lm_res)

    coeffs = lm_res['coefficients'].values()
    CI_m   = rpy.r.confint(lm_res['call'], 'x',           level=0.95)[0]
    CI_int = rpy.r.confint(lm_res['call'], '(Intercept)', level=0.95)[0]
    f.write('slope (95%% CI)\n%.4G'%(coeffs[0]) + '\n')
    f.write('[%.4G,%.4G]'%(CI_m[0],CI_m[1]) + '\n')

    f.write('intercept (95%% CI)\n%.4G'%(coeffs[1]) + '\n')
    f.write('[%.4G,%.4G]'%(CI_int[0],CI_int[1]) + '\n')

    #f.write('R-sq\n%.4G'%(lm_sum['r.squared']) + '\n')
    f.write('R-sq-adj\n%.4G'%(lm_sum['adj.r.squared']) + '\n')
    f.write('f-stat\n%.4G'%(lm_aov['F value'][0]) + '\n')
    f.write('Pr(>F)\n%.4G'%(lm_aov['Pr(>F)'][0]) + '\n')

    

if __name__ == "__main__":
    batchFile, dataDir, analysisMode = initialize()

    '''
    if analysisMode == 'runANDanalyze' or analysisMode == 'runANDstop':
       #see java -X | more for documentation
       runCmd = 'java -jar radicalization.jar -b ' + batchFile
       #runCmd = 'java -Xms500m -Xmx2G -jar radicalization.jar -b ' + batchFile
       print 'Running: ' + runCmd
       if os.system(runCmd):
           raise Exception, 'Batch execution failed.'
    '''
    if analysisMode == 'runANDanalyze' or analysisMode == 'analyzeOnly':
       print 'Starting data analysis on directory:'
       print dataDir
       report = reportStats(dataDir)
       displayStats(dataDir, report)


