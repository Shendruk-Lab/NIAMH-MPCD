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
import curvature_2D
###########################################################
### Calculates and plots end-to-end for 2D
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
myMap=ed.plasma
FS=20
test = False

#if test == True: #does for only one act/chunk rather than all
def myfunc():
    cwd = os.getcwd() #should be within a chunk
    time = []
    curvT = []
    pos = []
    curvL = []
    for i in range(1,21): #assuming 20 repeats
        temptime = []
        tempcurvT = []
        temppos = []
        tempcurvL = []
        path = str(cwd+"/result"+str(i))
        os.chdir(path)
        curvature_2D.func(path, 100, 100, 24)
        Lname = path+'/curvature-L.dat'
        tname = path+'/curvature-time.dat'
        Lfile = open(Lname,"r")
        tfile = open(tname,"r")
        for j in range(2): #toss headers
            Lline = Lfile.readline()
            tline = tfile.readline()
        for j in range(24):
            Lline = Lfile.readline()
            position, curv, std = Lline.split()
            temppos.append(float(position))
            tempcurvL.append(float(curv))
        pos = temppos #constant
        curvL.append(tempcurvL)
        for j in range(200):
            tline = tfile.readline()
            t, curv, std = tline.split()
            temptime.append(500*float(t))
            tempcurvT.append(float(curv))
        time = temptime
        curvT.append(tempcurvT)
    curvLarr = np.asarray(curvL)
    curvTarr = np.asarray(curvT)
    #average each
    curvLav = np.average(curvLarr, axis=0)
    curvTav = np.average(curvTarr, axis=0)    
    #get std and stderr for each
    Lstd = np.std(curvLarr, axis=0)
    Tstd = np.std(curvTarr, axis=0)
    Lerr = Lstd/np.sqrt(np.shape(curvLarr)[0])
    Terr = Tstd/np.sqrt(np.shape(curvTarr)[0])
    
    # write out 
    # Position one first
    bend,act,chunks = gyrationTensor_2D.read(cwd)
    name = str(cwd+"/curvL_a"+str(act)+"_c"+str(chunks)+".dat")
    outfile = open(name, 'w')
    outfile.write("Average curvature at each monomer position with 20 repeats\n")
    outfile.write("position\tcurv\tstd\tsterr\n")
    for t in range(len(pos)):
        outfile.write("%e\t%e\t%e\t%e\n"%(pos[t],curvLav[t],Lstd[t],Lerr[t]))	
    outfile.close()
    # Now time one
    name = str(cwd+"/curvT_a"+str(act)+"_c"+str(chunks)+".dat")
    outfile = open(name, 'w')
    outfile.write("Average curvature with time, also averaged over 20 repeats \n")
    outfile.write("time\tcurv\tstd\tsterr\n")
    for t in range(len(time)):
        outfile.write("%e\t%e\t%e\t%e\n"%(time[t],curvTav[t],Tstd[t],Terr[t]))	
    outfile.close()    
    # plot 2!!! figs, pos here
    fig,ax = plt.subplots(1)
    figdir = '/home/s1954660/Desktop/summer/analysis/curvPos'
    figname = str(figdir+"/curvL_b"+str(bend)+"_c"+str(chunks)+".pdf")
    ed.errorbar_fill(np.asarray(pos,dtype=float),curvLav,yerr=Lstd)
    title = str("Curvature at each monomer index. Bend: "+str(bend)+", chunks: "+str(chunks))
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(r"Monomer index",fontsize=FS)
    ax.set_ylabel(r"Average curvature",fontsize=FS)
    plt.savefig(figname,format='pdf', dpi = 'figure')
    plt.close(fig)
    # now time one
    fig,ax = plt.subplots(1)
    figdir = '/home/s1954660/Desktop/summer/analysis/curvTime'
    figname = str(figdir+"/curvT_b"+str(bend)+"_c"+str(chunks)+".pdf")
    ed.errorbar_fill(np.asarray(time,dtype=float),curvTav,yerr=Tstd)
    title = str("Average curvature over time. Bend: "+str(bend)+", chunks: "+str(chunks))
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(r"Time (simulation timesteps)",fontsize=FS)
    ax.set_ylabel(r"Average curvature",fontsize=FS)
    plt.savefig(figname,format='pdf', dpi = 'figure')
    plt.close(fig)

if test == True:
    myfunc()

else: #same but over all params (?!)
    # start in folder cluster01
    clustdir = os.getcwd()
    #datedirlist = os.listdir(clustdir)
    #for date in datedirlist: #go through date folders
    #    datedir = os.path.join(clustdir, date)
        #go through activities
    actdirlist = os.listdir(clustdir)
    for act in actdirlist:
        if 'act' in act:
            actdir = os.path.join(clustdir, act)
            #go through chunks
            chunkdirlist = os.listdir(actdir)
            for chunk in chunkdirlist:
                chunkdir = os.path.join(actdir, chunk)
                os.chdir(chunkdir)
                myfunc()