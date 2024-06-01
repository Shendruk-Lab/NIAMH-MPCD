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
mode = 1
radG = True

def myfunc():
    # os.chdir(beginning path if not already navigated here)
    time = []
    reelist = []
    cwd = os.getcwd()
    print(cwd)
    for i in range(1,21):
        temptime=[]
        tempree=[]
        path = str(cwd+"/result"+str(i))
        os.chdir(path)
        # create gyration.dat for that file
        gyrationTensor_2D.func(cwd)
        # access gyration.dat
        name = path+'/gyration.dat'
        file=open(name, "r")
        line=file.readline() #Gyration tensor measures and end to end distance with time 
        line=file.readline() #time	rg	rg_0	rg_1	k^2	ree
        for j in range(200):
            line=file.readline()
            t,rg,rg_0,rg_1,k2,ree = line.split()
            temptime.append(500*float(t))
            if radG==True:
                tempree.append(float(rg))
            else:
                tempree.append(float(ree))
        time=temptime
        #print(tempree)
        reelist.append(tempree)
    #print(time)
    #print(reelist)
    reearr = np.asarray(reelist,dtype=float)
    print(reearr)
    # then want to average along the correct!! axis
    print(np.shape(reearr))
    reeav = np.average(reearr, axis=0)
    reestd = np.std(reearr, axis=0)
    reeerr = np.std(reearr, axis=0)/np.sqrt(np.shape(reearr)[0])
    print("should be 20: ",np.shape(reearr)[0])
    print(reeerr)
    print(reeav)

    # write out everything in case need to come back to
    outfile = open(cwd+"/endToEnd.dat", 'w')
    outfile.write("Average end to end distance with time \n")
    outfile.write("time\tree\tstd\tsterr\n")
    for t in range(len(time)):
        outfile.write("%e\t%e\t%e\t%e\n"%(time[t],reeav[t],reestd[t],reeerr[t]))	
    outfile.close()

    fig,ax = plt.subplots(1)
    #plt.plot(time, reeav)
    figdir = '/home/s1954660/Desktop/for_tyler/rg/'
    bend,activity,chunks = gyrationTensor_2D.read(cwd)
    figname=str(figdir+"endToEnd_b"+str(activity)+"_c"+str(chunks)+".pdf")
    ed.errorbar_fill(np.asarray(time,dtype=float),reeav,yerr=reestd) #probs calculate an error
    title = str("Radius of gyration averaged over 20 runs, with standard devation. Act: "+str(activity)+". Chunks: "+str(chunks)+".")
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(r"Time (simulation timesteps)",fontsize=FS)
    ax.set_ylabel(r"Average radius of gyration (simulation cells)",fontsize=FS)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    plt.savefig(figname,format='pdf', dpi = 'figure')
    plt.close(fig)
    
if mode == 0:
    myfunc()

elif mode ==1: #same but over all params (?!)
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
                
else: #for comparison/ multiple lines on the one graph
    print("yes")
    ext = os.path.join(os.getcwd(), 'act-4','chunks-1') #cluster05
    # os.chdir(beginning path if not already navigated here)
    con = os.path.join(os.getcwd(), 'act--4','chunks-1')
    zero = os.path.join(os.getcwd(), 'act-0')
    # have things i can read in
    names = ['act-4/chunks-2','act-4/chunks-3','act--4/chunks-3','act-4/chunks-4']
    labels = ['diblock copolymer', 'triblock (ext. outer)', 'triblock (con. outer)', 'tetrablock']#, 'passive polymer']
    bigtime = []
    bigrg = []
    bigstd = []
    def add(fname,i):
        time = []
        radg = []   
        weestd = []
        rname = os.path.join(os.getcwd(), names[i],'endToEnd.dat')
        file = open(rname,'r')
        #toss headers
        for j in range(2):
            line=file.readline()
        while file:
            line = file.readline()
            if (not line):
                break
            line=line.split()
            t = line[0]
            #if float(t)>30000:
            #    break
            e = line[1]
            #if i<2:
            std = line[2]
            time.append(float(t))
            #else:
            #    time.append(25*float(t))
            radg.append(float(e))
            #if i<2:
            weestd.append(float(std))
        file.close()
        #print(time)
        bigtime.append(time)
        bigrg.append(radg)
        #if i<2:
        bigstd.append(weestd)
        #else:
        #    bigstd.append(0.0)
    for i in range(len(names)):
        fname = os.path.join(names[i], 'endToEnd.dat')
        add(fname,i)
    lines = []
    fig,ax = plt.subplots(1)
    for i in range(len(names)):
        lines.append(ed.errorbar_fill(np.asarray(bigtime[i],dtype=float),np.asarray(bigrg[i]),yerr=np.asarray(bigstd[i])))
    ax.legend(labels, fontsize=12, loc=2)
    ax.set_xlabel(r"time, $t_{\mbox{mpcd}}$",fontsize=FS)
    ax.set_ylabel(r"Radius of gyration",fontsize=FS)
    ax.set_title(r"Comparison of different copolymers, activity = 4",fontsize=FS)
    plt.savefig('/home/s1954660/Desktop/for_tyler/rg'+"/ChunksComparisonA4.pdf",format='pdf', dpi = 'figure')
    plt.savefig('/home/s1954660/Desktop/for_tyler/rg'+"/ChunksComparisonA4.png",format='png', dpi = 'figure')
    plt.close(fig)