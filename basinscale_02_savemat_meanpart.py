from __future__ import print_function

import os.path as p
import subprocess as sub
from builtins import range

import glob

import numpy as np

import scipy.io as scio

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))
root_path = p.dirname(self_path)
'''
import sys
sys.path.append(p.join(root_path,'striations_3layers/try06_S1'))

print('self_path=',self_path)
print('root_path=',root_path)
'''
##### plot wind forcing
# basinscale_3layer_s0
xlen0 = [2880, 3840, 5760, 7680];
ylen0 = [3840, 3840, 3840, 3840];
nx0 =   [ 97, 129, 193, 257];
ny0 =   [129, 129, 129, 129];

# basinscale_3layer_re15km
xlen1 = [ 960, 1920, 2880, 3840, 4800, 5760];
ylen1 = [3840, 3840, 3840, 3840, 3840, 3840];
nx1 = [ 64, 128, 192, 256, 320, 384];
ny1 = [256, 256, 256, 256, 256, 256];

# basinscale_3layer_inceaseRE
xlen2 = [3840, 3840, 3840, 3840];
ylen2 = [3840, 3840, 3840, 3840];
nx2 = [129, 192, 256, 384];
ny2 = [129, 192, 256, 384];

layers = 3
#grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

for s in range(1,7):
    print(s)
    xlen = xlen1[s-1]*1.e3
    ylen = ylen1[s-1]*1.e3
    nx = nx1[s-1]
    ny = ny1[s-1]
    grid  = aro.Grid(nx, ny, layers, xlen/nx, ylen/ny)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ##### read
    v_files =   sorted(glob.glob('./s1' + str(s) + "/output/av.v.*"))
    u_files =   sorted(glob.glob('./s1' + str(s) + "/output/av.u.*"))
    eta_files = sorted(glob.glob('./s1' + str(s) + "/output/av.eta.*"))
    h_files =   sorted(glob.glob('./s1' + str(s) + "/output/av.h.*"))

    print( "len(v_files) = ", len(v_files))
    print( "len(u_files) = ", len(u_files))

    y1 = len(v_files) - 5
    y2 = len(v_files) 
  
    ##### plot each state of the run
    v1 = np.zeros((layers, nx,ny+1))
    u1 = np.zeros((layers, nx+1,ny))
    eta1 = np.zeros((1, nx,ny))
    h1 = np.zeros((layers, nx,ny))
    #print(type(u1))
    for i in range(y1,y2,1):
        print(i)
        v = aro.interpret_raw_file(v_files[i], grid.nx, grid.ny, grid.layers)
        u = aro.interpret_raw_file(u_files[i], grid.nx, grid.ny, grid.layers)
        eta = aro.interpret_raw_file(eta_files[i], grid.nx, grid.ny, grid.layers)
        h = aro.interpret_raw_file(h_files[i], grid.nx, grid.ny, grid.layers)
    
        #print(type(u))
        print(eta.shape) 
        print(h.shape) 
        #print(u)
        v[np.isnan(v)] = 0
        u[np.isnan(u)] = 0
        eta[np.isnan(eta)] = 0
        h[np.isnan(h)] = 0
        #print(u)
    
        v1 += v.transpose(0,2,1)
        u1 += u.transpose(0,2,1)
        eta1 += eta.transpose(0,2,1)
        h1 += h.transpose(0,2,1)

    #print(u1)
    
    dataNew1 = './s1' + str(s) + '/output_v_mean_141-150.mat'
    dataNew2 = './s1' + str(s) + '/output_u_mean_141-150.mat'
    dataNew3 = './s1' + str(s) + '/output_eta_mean_141-150.mat'
    dataNew4 = './s1' + str(s) + '/output_h_mean_141-150.mat'

    scio.savemat(dataNew1, {'v': v1/10})
    scio.savemat(dataNew2, {'u': u1/10})
    scio.savemat(dataNew3, {'eta': eta1/10})
    scio.savemat(dataNew4, {'h': h1/10})

