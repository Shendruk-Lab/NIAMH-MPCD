import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker
from pylab import *
from numpy import ma
from subprocess import call
from numpy import linalg as LA
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import os
import csv
import sys
import numpy as np
import math
import gyrationTensor_2D
###########################################################
### Calculates and plots end-to-end for 2D
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
myMap=ed.plasma
FS=15
saveloc = '/home/s1954660/Desktop/summer/analysis/velocities'

# find cm using gyrationTensor_2D.py
# don't necessarily want to be averaging any velocities - instead
# can plot EVERY curv and EVERY velocity for each repeat.
# HOWEVER I can see this being super noisy
# all on the same graph for each chunk type? probs not?


def myfunc(dir, prevwd):
    # i think i am going to assume that the unwrap thing is for periodic boundaries
    cm, dt = gyrationTensor_2D.cmfunc(dir, prevwd)
    sp = np.shape(cm) #(2, 80020)
    vel = np.zeros_like(cm)
    print(np.shape(vel))
    #print(np.shape(cm))
    # forward difference for t=0
    # backward difference for t=-1 
    for d in range(2):
        vel[d][0]=(cm[d][1]-cm[d][0])/dt
        vel[d][-1]=(cm[d][-1]-cm[d][-2])/dt
    # centered difference for the rest
    for t in range(1,sp[-1]-1):
        for d in range(2):
            vel[d][t]=(cm[d][t+1]-cm[d][t-1])/(2*dt)
            
    # just get magnitude? use the norm thingy
    speed = np.linalg.norm(vel, axis=0)
    print(np.shape(speed))
    for i in range(sp[-1]):
        allvels.append(speed[i])
    # now read in rg and curv dat files.
    rgname = 'gyration.dat'
    curvname = 'curvature-time.dat'
    rgs = []
    fname = os.path.join(os.getcwd(), rgname)
    f = open(fname,'r')
    #toss headers
    line=f.readline()
    line=f.readline()
    while f:
        line = f.readline()
        if (not line):
            break
        line = line.split()
        allrgs.append(float(line[1]))
    f.close()
    curvs = []    
    fname = os.path.join(os.getcwd(), curvname)
    f = open(fname,'r')
    #toss headers
    line=f.readline()
    line=f.readline()
    while f:
        line = f.readline()
        if (not line):
            break
        line = line.split()
        allcurvs.append(float(line[1]))
    f.close()
    
def myplot(allvels,allcurvs,allrgs,dir):
    bend,act,chunks = gyrationTensor_2D.read(dir)
    rgname = str('velVsRg_a'+str(act)+'_c'+str(chunks)+'.pdf')
    curvname = str('velVsCurv_a'+str(act)+'_c'+str(chunks)+'.pdf')
    #rg plot
    fig,ax = plt.subplots(1)
    figname = os.path.join(saveloc, rgname)
    plt.plot(allrgs, allvels,'o', markersize=1)
    title = str("Velocity of polymer against radius of gyration. Act: "+str(act)+", chunks: "+str(chunks))
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(r"Radius of gyration",fontsize=FS)
    ax.set_ylabel(r"Velocity of polymer",fontsize=FS)
    plt.savefig(figname,format='pdf', dpi = 'figure')
    plt.close(fig)
    #curv plot
    fig,ax = plt.subplots(1)
    figname = os.path.join(saveloc, curvname)
    plt.plot(allcurvs, allvels,'o', markersize=1)
    title = str("Velocity of polymer against curvature. Act: "+str(act)+", chunks: "+str(chunks))
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(r"Average curvature along backbone",fontsize=FS)
    ax.set_ylabel(r"Velocity of polymer",fontsize=FS)
    plt.savefig(figname,format='pdf', dpi = 'figure')
    plt.close(fig)
    
    


#start in cluster05, go through each activity
clustdir = os.getcwd()
actdirlist = os.listdir(clustdir)
for act in actdirlist:
    if 'act' in act:
        actdir = os.path.join(clustdir, act)
        #go through chunks
        chunkdirlist = os.listdir(actdir)
        for chunk in chunkdirlist:
            chunkdir = os.path.join(actdir, chunk)
            resdirlist = os.listdir(chunkdir)
            allvels = []
            allcurvs = []
            allrgs = []
            for res in resdirlist:
                if 'result' in res:
                    resdir = os.path.join(chunkdir, res)
                    os.chdir(resdir)
                    myfunc(os.getcwd(), chunkdir)
            myplot(allvels,allcurvs,allrgs,chunkdir)