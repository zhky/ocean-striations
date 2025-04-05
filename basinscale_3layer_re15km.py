from __future__ import print_function

import os.path as p
import subprocess as sub
from builtins import range

import glob

import numpy as np
import random

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)

print('self_path=',self_path)
print('root_path=',root_path)


# The examples
#
# N.B. These configurations are equivalent to several of those in the test
# suite, but run at a higher resolution.

def wind_asymmetric(X, Y): #the wind input file is in Nm^-2 
    """ follow shevchenko and Berloff, 2015 """
    wind_forcing = np.zeros(Y.shape)

    x = X[0,::]
    y = Y[::,0]
    #print(x)
    #print(x)

    for i in range(1,len(x)):
        y0 = 0.4*Y.max() + 0.2*x[i-1]
        #print(y-y0)

        for j in range(1,len(y)):
            #print(y[j-1]-y0)
            if y[j-1] < y0:
                wind_forcing[j-1,i-1] = -1.8*np.pi*0.03*(np.sin(np.pi*y[j-1]/(y0)))
            else:
                wind_forcing[j-1,i-1] = 2.22*np.pi*0.03*(np.sin(np.pi*(y[j-1]-y0)/(Y.max()-y0)))

    plt.figure()
    plt.plot(wind_forcing[:,1], Y[:,1]/1e3)
    plt.ylabel('Latitude (km)')
    plt.xlabel('Wind forcing (N/m^2)')
    plt.savefig('wind_forcing.png',bbox_inches='tight',dpi=200)
    plt.close()

    return wind_forcing

def sponge_h1(X, Y):
    """Produce the sponge file of layer 1"""

    sponge_h = 600*np.ones(X.shape, dtype=np.float64)

    plt.figure()
    plt.pcolormesh(X,Y,sponge_h,shading='auto')
    plt.colorbar()
    plt.savefig('sponge_h.png', dpi=150)
    plt.close()

    return sponge_h

def sponge_h2(X, Y):
    """Produce the sponge file of layer 2"""
    sponge_h = 1400*np.ones(X.shape, dtype=np.float64)

    return sponge_h

def sponge_h3(X, Y):
    """Produce the sponge file of layer 3"""
    sponge_h = 2000*np.ones(X.shape, dtype=np.float64)

    return sponge_h


def sponge_h_timescale(X, Y):
    """Produce the sponge timescale file"""
    sponge_h_timescale = np.zeros(X.shape, dtype=np.float64)
    sponge_h_timescale[Y<5*Y.min()] = 1/(2.*30.*86400.) # 2 month relaxation time
    sponge_h_timescale[Y>Y.max()-5*Y.min() ] = 1/(2.*30.*86400.)

    plt.figure()
    plt.pcolormesh(X,Y,sponge_h_timescale*86400.*30.,shading='auto')
    plt.colorbar()
    plt.savefig('sponge_h_timescale.png',dpi=150)
    plt.close()

    return sponge_h_timescale

def wetmask(X, Y):
    """The wet mask."""

    # start with land everywhere and carve out space for water
    wetmask = np.ones(X.shape, dtype=np.float64)

    # clean up the edges
    wetmask[ 0, :] = 0
    wetmask[-1, :] = 0
    wetmask[ :, 0] = 0
    wetmask[ :,-1] = 0

    plt.figure()
    plt.pcolormesh(X/1e3, Y/1e3, wetmask, shading='auto', cmap='Greys_r')
    plt.colorbar()
    #plt.xlim(0,1500)
    #plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')
    plt.savefig('wetmask.png', dpi=150, bbox_inches='tight')
    plt.close()

    return wetmask

def coriolis(X,Y):
    r = Y - Y.max()/2
    f0 = 7.3e-5 # at middle latitude 
    beta = f0*np.cos(np.pi*40/180)/6371e3

    f = f0 + beta*r

    plt.pcolormesh(X, Y, f, shading='auto')
    plt.colorbar()
    plt.savefig('coriolis.png')
    plt.close()
    return f

def wind_asymmetric_steady(xlen, ylen, nx, ny, nTimeSteps, dt, botDrag, simulation=None):
    layers = 3
    dumpFreq = 360*24*3600
    avFreq = 360*24*3600*2

    with working_directory(p.join(self_path, 
        "{0}".format(simulation))):
        sub.check_call(["rm","-rf","output/","input/","*.mat"])
        drv.simulate(initHfile=[600., 1400., 3000.],
                zonalWindFile=[wind_asymmetric],
                wind_mag_time_series_file = [0.1],
                wetMaskFile = [wetmask],
                fUfile=[coriolis],
                fVfile=[coriolis],
                botDrag = botDrag,
                RelativeWind = 'no',
                exe='aronnax_external_solver',
                nx=nx, ny=ny, layers=layers,
                dx=xlen/nx, dy=ylen/ny,
                dt = dt,
                nTimeSteps = nTimeSteps,  #  total years (nTimeSteps*dt)
                dumpFreq = dumpFreq,   # yearly snapshot output
                avFreq = avFreq     # yearly average output
                )

    
if __name__ == '__main__':


    # increase horizontal resolution
    #wind_asymmetric_steady( 960e3,3840e3, 64,256,150*360*24*4,900,4.e-8,'s11')
    #wind_asymmetric_steady(1920e3,3840e3,128,256,150*360*24*6,600,4.e-8,'s12')
    #wind_asymmetric_steady(2880e3,3840e3,192,256,150*360*24*6,600,4.e-8,'s13')
    #wind_asymmetric_steady(3840e3,3840e3,256,256,150*360*24*6,600,4.e-8,'s14')
    #wind_asymmetric_steady(4800e3,3840e3,320,256,150*360*24*6,600,4.e-8,'s15')
    wind_asymmetric_steady(5760e3,3840e3,384,256,150*360*24*6,600,4.e-8,'s16')
 
