from __future__ import print_function

import os.path as p
import subprocess as sub
from builtins import range

import glob

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import aronnax as aro
import aronnax.driver as drv
from aronnax.utils import working_directory

self_path = p.dirname(p.abspath(__file__))


##### initial parameters

# basinscale_3layer_re15km
xlen1 = [ 960, 1920, 2880, 3840, 4800];
ylen1 = [3840, 3840, 3840, 3840, 3840];
nx1 = [ 64, 128, 192, 256, 320];
ny1 = [256, 256, 256, 256, 256];

layers = 3

for s in range(4,5):
    
    xlen = xlen1[s-1]*1.e3
    ylen = ylen1[s-1]*1.e3
    nx = nx1[s-1]
    ny = ny1[s-1]

    grid = aro.Grid(nx, ny, layers, xlen / nx, ylen / ny)

    ##### plot the mean state

    v_files = sorted(glob.glob('./s1'+str(s)+'/output/snap.v.*'))
    u_files = sorted(glob.glob('./s1'+str(s)+'/output/snap.u.*'))

    
    print( "len(v_files) = ", len(v_files))
    print( "len(u_files) = ", len(u_files))
    
    y = 149 # output is yearly
    
    v = aro.interpret_raw_file(v_files[y], grid.nx, grid.ny, grid.layers)
    u = aro.interpret_raw_file(u_files[y], grid.nx, grid.ny, grid.layers)
    #print("u.shape = ", np.shape(u))
    #print("v.shape = ", np.shape(v))
    #print(u)
    #print(v)
    
    
    # plot ---------------------------------------------
    fig = plt.figure(figsize=(14,4),facecolor='white')
    lay = (0,1,2)

    for i in lay:
    
        ax = fig.add_subplot(1,3,i+1)
        X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)
        plt.streamplot(X,Y,u[i,:,:-1],v[i,:-1,:],color='k',density=1.5)
                 
        X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)      
        colour_lim = u.max()*100*0.6
        im = plt.pcolormesh(X,Y,u[i,:,:-1]*100., cmap='RdBu_r', shading='auto'
            ,vmin = -colour_lim, vmax = colour_lim)
        im.set_edgecolor('face')
        
        CB = plt.colorbar()
        CB.set_label('x component of velocity (cm / s)')

        #plt.axes().set_aspect('equal')
        plt.xlabel('x coordinate (km)')
        plt.ylabel('y coordinate (km)')
        plt.title( str(y) + 'th year, layer' + str(lay[i]+1) )
        plt.tight_layout()
        
    plt.show()  
        #plt.savefig('./try06/layer_average.png'.format(u_files[i][-10:]), dpi=150,
        #    bbox_inches='tight')
    #plt.close()


############################################################################
'''    
    
##### plot each state of the run
v_files = sorted(glob.glob("./try06/output/snap.v.*"))
u_files = sorted(glob.glob("./try06/output/snap.u.*"))

print( "len(v_files) = ", len(v_files))
print( "len(u_files) = ", len(u_files))

for i in range(len(v_files)):
    v = aro.interpret_raw_file(v_files[i], grid.nx, grid.ny, grid.layers)
    u = aro.interpret_raw_file(u_files[i], grid.nx, grid.ny, grid.layers)

    # ---------------------------------------------
    plt.figure()
    ln = 0   # layer1
    X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)
    plt.streamplot(X,Y,u[ln,:,:-1],v[ln,:-1,:],color='k',density=1)
         
    X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)        
    colour_lim = 1
    im = plt.pcolormesh(X,Y,u[ln,:,:-1]*100., cmap='RdBu_r'
        ,vmin = -colour_lim, vmax = colour_lim)
    im.set_edgecolor('face')
    CB = plt.colorbar()
    CB.set_label('x component of velocity (cm / s)')
    plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')
    plt.title('layer 1: timestep {0}'.format(u_files[i][-10:]))
        
    plt.savefig('./try06/figure/layer1_{0}.png'.format(u_files[i][-10:]), dpi=150,
        bbox_inches='tight')
       
    plt.close()


    # ----------------------------------------------
    plt.figure()
    ln = 1  # layer2
    X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)
    plt.streamplot(X,Y,u[ln,:,:-1],v[ln,:-1,:],color='k',density=1)
         
    X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)   # what is xp1, yp1        
    colour_lim = .5
    im = plt.pcolormesh(X,Y,u[ln,:,:-1]*100., cmap='RdBu_r'
        ,vmin = -colour_lim, vmax = colour_lim)
    im.set_edgecolor('face')
    CB = plt.colorbar()
    CB.set_label('x component of velocity (cm / s)')
    plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')
    plt.title('layer 2: timestep {0}'.format(u_files[i][-10:]))
        
    plt.savefig('./try06/figure/layer2_{0}.png'.format(u_files[i][-10:]), dpi=150,
        bbox_inches='tight')
        
    plt.close()
    
    # ----------------------------------------------
    plt.figure()
    ln = 2  # layer3
    X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)
    plt.streamplot(X,Y,u[ln,:,:-1],v[ln,:-1,:],color='k',density=1)
         
    X,Y = np.meshgrid(grid.x/1e3, grid.y/1e3)   # what is xp1, yp1        
    colour_lim = .1
    im = plt.pcolormesh(X,Y,u[ln,:,:-1]*100., cmap='RdBu_r'
        ,vmin = -colour_lim, vmax = colour_lim)
    im.set_edgecolor('face')
    CB = plt.colorbar()
    CB.set_label('x component of velocity (cm / s)')
    plt.axes().set_aspect('equal')
    plt.xlabel('x coordinate (km)')
    plt.ylabel('y coordinate (km)')
    plt.title('layer 3: timestep {0}'.format(u_files[i][-10:]))
        
    plt.savefig('./try06/figure/layer3_{0}.png'.format(u_files[i][-10:]), dpi=150,
        bbox_inches='tight')
        
    plt.close()


try:
    sub.check_call(["convert", "-delay", "30", "-loop", "0", "./try06/figure/layer1_*.png", "./try06/figure/{0}.gif".format('layer1')])
except:
    print("failed to make animation")
    
try:
    sub.check_call(["convert", "-delay", "30", "-loop", "0", "./try06/figure/layer2_*.png", "./try06/figure/{0}.gif".format('layer2')])
except:
    print("failed to make animation")
        

try:
    sub.check_call(["convert", "-delay", "30", "-loop", "0", "./try06/figure/layer3_*.png", "./try06/figure/{0}.gif".format('layer3')])
except:
    print("failed to make animation")

'''
