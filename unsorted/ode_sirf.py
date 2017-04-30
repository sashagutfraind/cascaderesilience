'''
    Plotting of the SIRF model
'''

import matplotlib
import sys
#for running scripts silently:
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('PS')
#if 'matplotlib.backends' not in sys.modules: matplotlib.use('pdf')
#import matplotlib.pylab as pylab
import pylab

import inspect
import numpy as np
import numpy.random as npr
import csv, time, os
import pdb
import scipy
import scipy.integrate

SAVEFIGS = True

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())

'''
class reportFile:
    f = None
    def __init__(self,fname=None):
        if fname == None:
            fname = 'output/report_'+timeNow()+'.csv'
        self.f = open(fname, 'wb')
    def close(self):
        if self.f!= None:
            self.f.close()
    def writeln(self,s=''):
        self.f.write(str(s) + os.linesep)
        print s

reportF = reportFile()
'''

def contour_plot(Xs, Ys, Xname, Yname, odeParams, initialVals, caption='R infinity', filename=None):
    #Xmesh, Ymesh = np.meshgrid(Xs, Ys) #for some reason, this doesn't work if Xs and Ys have different lengths!

    Zmesh = np.zeros(shape=(len(Xs), len(Ys)), dtype=np.double)
    for xindex, x in enumerate(Xs):
        for yindex, y in enumerate(Ys):
            params = dict(odeParams)
            params[Xname] = x
            params[Yname] = y

            final_epidemic = final_values(odeParams=params.items(), initialVals=initialVals)
    
            Zmesh[xindex, yindex] = final_epidemic['R']

    #strangely enough, mpl uses Ys for the first index!
    Zmesh = Zmesh.transpose()
    
    import matplotlib.pyplot as plt
    plt.figure()
    #CS = plt.contour(Xs, Ys, Zmesh)
    CS = plt.contour(Xs, Ys, Zmesh, levels=np.arange(0,1.0, 0.05))
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title(caption)
    plt.xlabel(Xname, fontsize=20)
    plt.ylabel(Yname, fontsize=20)

    if filename==None:
        filename = 'output/contour_'+timeNow()+'.pdf'
    plt.savefig(filename)
    pylab.show()
    #see more features in:
    #http://matplotlib.sourceforge.net/plot_directive/mpl_examples/pylab_examples/contour_demo.py



def dxdt(x, t, params):
#the code with a 
    S = x[0]
    I = x[1]
    R = x[2]
    F = x[3]
    N = S + I + R + F

    for pname,pval in params:
        exec pname + '=%.20f'%pval  #this might lose precision
        #params must be a tuple
        #so they are given as a tuple of tuples (('param1',1.0), ('param2',22.))

    #pdb.set_trace()
    #fear_factor =  gamma * ((I/(1.*N))**(1+p))
    fear_factor =  I > p and gamma or 0.0000001
    
    rateS = -beta*S*(I/(1.*N)) - S*fear_factor
    rateI = +beta*S*(I/(1.*N)) - rho * I
    rateR =                    + rho * I
    rateF =                    + S*fear_factor

    return [rateS, rateI, rateR, rateF]


def final_values(odeParams, initialVals, systemFunc=dxdt, tmax=50., dt=.001):
    x0          = [val for varName,val in initialVals]
    times       = scipy.arange(0., tmax, dt, dtype=scipy.float32)
    xOrbit      = scipy.integrate.odeint(systemFunc, x0, times, args=(odeParams,))
  
    ret = {}
    for varNum,varInfo in enumerate(initialVals):
        varName=varInfo[0]
        ret[varName] = xOrbit[-1,varNum]

    return ret


def plot_equilibrium_pts():
    import mpl_toolkits.mplot3d.axes3d as a3d
    import matplotlib.pyplot as plt

    RGs = np.arange(0, 2.5, 0.05)
    RIs = np.arange(0, 2.5, 0.05)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #X, Y, Z = a3d.get_test_data(0.05)
    
    RGbound = lambda rg,ri: 1./(1+rg)
    RIbound = lambda rg,ri: ri/(1+ri)

    #lowerZ =  

    ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

    plt.show()

    pylab.savefig('output/equils_'+timeNow(), dpi=300)

def plot_equilibrium_pts2d():
    RGs = np.arange(1, 4.5, 0.1)
    RIs = np.arange(1, 4.5, 0.1)

    RGbound = lambda rg: 1./(1+rg)
    RIbound = lambda ri: ri/(1+ri)

    lowerZ = map(RGbound, RGs)
    upperZ = map(RIbound, RIs)

    pylab.figure()
    pylab.hold(True)
    pylab.plot(RGs, lowerZ, 'b-^')
    pylab.plot(RIs, upperZ, 'r-')
    pylab.xlim(min(RIs), max(RIs)*0.99)
    pylab.ylim(0,1.0)

    pylab.xlabel(r'$r_G = r_I$', fontsize=20)
    pylab.ylabel(r'$S_0$', fontsize=20)

    pylab.hold(False)

    pylab.text(x=1.3, y=0.9, s='government victory', fontsize=20)
    pylab.text(x=2.5, y=0.5, s='stalemate', fontsize=20)
    pylab.text(x=1.3, y=0.1, s='insurgent victory', fontsize=20)

    pylab.savefig('output/equils2d_'+timeNow(), dpi=300)


def plot_vs_t(odeParams, initialVals, systemFunc=dxdt, tmax=50., dt=.001, plottedSeries = None):
    #if plottedSeries == None:
    #    plottedSeries = [varName for varName,val in initialVals]
    x0          = [val for varName,val in initialVals]
    times       = scipy.arange(0., tmax, dt, dtype=scipy.float32)
    xOrbit      = scipy.integrate.odeint(systemFunc, x0, times, args=(odeParams,) )
    #xOrbitM     = {}
    #(xOrbitM[0], xOrbitM[1]) = purifyData(xOrbit,rescale=False)
    #xOrbitM[2] = getCapacity(xOrbitM[0], xOrbitM[1], m_)
   
    fig=pylab.figure()
    pylab.hold(True)
    curvetypes = ['r-', 'g-', 'b-', 'k-^']
    linewidths = [4.0, 2.5, 1.0, 0.5]
    for varNum,varInfo in enumerate(initialVals):
        varName=varInfo[0]
        #if varName not in plottedSeries:
        #    continue
        pylab.plot(times, xOrbit[:,varNum], curvetypes[varNum], label=varName, linewidth=linewidths[varNum])
    pylab.legend(loc='best', shadow=True)
    
    hregion     = tmax
    fregion     = np.max(xOrbit)
    pylab.axis([0, hregion, 0, fregion])
    paramstr = str(['%s=%.3f'%(varName,val) for varName,val in odeParams])
    #pylab.text(hregion/4., -fregion/10., paramstr)
    pylab.text(0.05*hregion, 0.7*fregion, paramstr, fontsize=14)
    
    xOstring = 'Initial Size:\n'+str(['%s=%.3f'%(varName,val) for varName,val in initialVals])
    pylab.text(0.05*hregion, 0.8*fregion, xOstring, fontsize=14)
    pylab.xlabel('time t', position=(1.,1.), fontsize=14)
    if SAVEFIGS:
        pylab.savefig('output/vs_T'+timeNow(), dpi=300)

    pylab.hold(False)
    pylab.show()

def purifyData(xOrbit, rescale=True):
#rescale and remove unphysical negative values
    index  = 0
    maxidx = scipy.shape(xOrbit)[1]
    while index < maxidx:
        orbit = xOrbit[:,index]
        if rescale:
            orbit = orbit/max(max(orbit),0.00001)
        for jndex, pt in enumerate(orbit):
            orbit[jndex] = max(0., pt)
        xOrbit[:,index] = orbit
        index += 1
    return (xOrbit[:,0], xOrbit[:,1])


if __name__ == '__main__':
    '''
    plot_vs_t(initialVals=(('S',0.9999),('I',0.0001),('R',0.),('F',0.)),
              odeParams=(('beta',0.1),('rho',1),('gamma',30),('p', 0.1),), 
              )

    plot_vs_t(initialVals=(('S',0.9999),('I',0.0001),('R',0.),('F',0.)),
              odeParams=(('beta',1.2),('rho',1),('gamma',30),('p', 0.1),), 
              )

    plot_vs_t(initialVals=(('S',0.9999),('I',0.0001),('R',0.),('F',0.)),
              odeParams=(('beta',2.5),('rho',1),('gamma',30),('p', 0.1),), 
              )

    '''
    contour_plot(
            Xs = np.arange(0., 5., 0.2),
            #Ys = np.arange(0., 50, 1.0),
            #Xs = np.arange(0., 100., 1.),
            Ys = np.arange(0., 30,  2.),
            Xname = 'beta',
            Yname = 'gamma',
            #odeParams = (('beta',2.),('rho',1),('gamma',1000),('p', 1.)),
            odeParams = (('beta',2.),('rho',1),('gamma',1000),('p', 0.1)),
            initialVals = (('S',0.9999),('I',0.0001),('R',0.),('F',0.)),
              )
