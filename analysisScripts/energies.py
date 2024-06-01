import matplotlib
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
import shendrukGroupFormat as ed
plt.style.use('shendrukGroupStyle')

cwd = os.getcwd() #should be energy_checks
FS=20

def add(fname):
    time = []
    energy = []   
    rname = os.path.join(fname, 'energy.dat') 
    file = open(rname,'r')
    #toss headers
    for j in range(13):
        line=file.readline()
    while file:
        line = file.readline()
        if (not line):
            break
        line=line.split()
        t = line[0]
        if float(t)>30000:
            break
        e = line[-1]
        time.append(float(t))
        energy.append(abs(float(e)-1))
    file.close()
    #print(time)
    bigtime.append(time)
    bigenergy.append(energy)
    #print(energy)
    #ax.plot(np.array(time), np.array(energy))
def avadd(fname): #for the ones with repeats ALSO ADD STDEV TO THIS 
    time = []
    energy = []    
    for i in range(5): #5 repeats
        tempt = []
        tempe = []
        rname = os.path.join(fname, str('result'+str(i+1)),'energy.dat')
        print(rname)
        file = open(rname,'r')
        #toss headers
        for j in range(13):
            line=file.readline()
        while file:
            line = file.readline()
            if (not line):
                break
            line=line.split()
            t = line[0]
            if float(t)>30000:
                break
            e = line[-1]
            tempt.append(float(t))
            tempe.append(abs(float(e)-1))
        file.close()
        time=tempt
        energy.append(tempe)
    stden = np.std(np.array(energy), axis=0)
    energy = np.average(np.array(energy), axis=0)
    #print(time)
    bigtime.append(time)
    enlist=[energy[i] for i in range(np.shape(energy)[0])]
    bigenergy.append(enlist)
    bigstd.append(stden)

    
#names = ['act-0', 'act-0.1', 'act-0.2', 'act-0.4', 'act-1', 'act-2', 'act--0.1', 'act--0.2', 'act--0.4', 'act--1', 'act--2']
#labels = ['strength 0', 'strength 0.1', 'strength 0.2', 'strength 0.4', 'strength 1', 'strength 2', 'strength -0.1', 'strength -0.2', 'strength -0.4', 'strength -1', 'strength -2']
names = ['x0.5/a1_sysx0.5','x0.5/zM-1000','x0.5/zM-5000','x1/a1_sysx1','x1/zM-1000','x1/zM-5000']#,'ts-1','ts-2','ts-3'
labels = ['nothing extra, 50x50', 'zero mom f1000, 50x50','zero mom f5000, 50x50','nothing extra, 100x100', 'zero mom f1000, 100x100','zero mom f5000, 100x100']# 'thermostat 1','thermostat 2','thermostat 3',
bigtime = []
bigenergy = []
bigstd=[0]
fig,ax = plt.subplots(1)
fname = os.path.join(cwd,names[0])
add(fname)
for i in range(1,len(names)): #11
    #print(cwd)
    fname = os.path.join(cwd,names[i])#,'energy.dat')
    #print(fname)
    #fig,ax = plt.subplots(1)
    #ed.errorbar_fill(np.array(time),np.array(avK),yerr=np.array(stdK))
	# ax.plot(time, array)
    if i==3:
        add(fname)
    else:
        avadd(fname)
#ax.plot(np.array(bigtime[0]),np.array(bigenergy[0]),label=labels[0], np.array(bigtime[1]),np.array(bigenergy[1]), np.array(bigtime[2]),np.array(bigenergy[2]), np.array(bigtime[3]),np.array(bigenergy[3]), np.array(bigtime[4]),np.array(bigenergy[4]), np.array(bigtime[5]),np.array(bigenergy[5]), np.array(bigtime[6]),np.array(bigenergy[6]))
# do the log bit
print(bigenergy)
lines = []
for i in range(len(names)):
    lines.append(ax.semilogy(np.array(bigtime[i]),np.array(bigenergy[i]),label=labels[i], lw=0.5))
#line0=ax.semilogy(np.array(bigtime[0]),np.array(bigenergy[0]),label=labels[0], lw=0.5)
#line1=ax.semilogy(np.array(bigtime[1]),np.array(bigenergy[1]),label=labels[1], lw=0.5)
#line2=ax.semilogy(np.array(bigtime[2]),np.array(bigenergy[2]),label=labels[2], lw=0.5)
print(lines[0])
#line3=ax.semilogy(np.array(bigtime[2]),np.array(bigenergy[2]),label=labels[2], lw=0.5)
#print(line3)
#line0, = ax.semilogy(np.array(bigtime[0]),np.array(bigenergy[0]),label=labels[0], lw=0.5)
#line1, = ax.semilogy(np.array(bigtime[1]),np.array(bigenergy[1]),label=labels[1], lw=0.5)
#line2, = ax.semilogy(np.array(bigtime[2]),np.array(bigenergy[2]),label=labels[2], lw=0.5)
#line3, = ax.semilogy(np.array(bigtime[3]),np.array(bigenergy[3]),label=labels[3], lw=0.5)
#line4, = ax.semilogy(np.array(bigtime[4]),np.array(bigenergy[4]),label=labels[4], lw=0.5)
#line5, = ax.semilogy(np.array(bigtime[5]),np.array(bigenergy[5]),label=labels[5], lw=0.5)
#line6, = ax.semilogy(np.array(bigtime[6]),np.array(bigenergy[6]),label=labels[6], lw=0.5)
#line7, = ax.semilogy(np.array(bigtime[7]),np.array(bigenergy[7]),label=labels[7], lw=0.5)
#line8, = ax.semilogy(np.array(bigtime[8]),np.array(bigenergy[8]),label=labels[8], lw=0.5)
#line9, = ax.semilogy(np.array(bigtime[9]),np.array(bigenergy[9]),label=labels[9], lw=0.5)
#line10, = ax.semilogy(np.array(bigtime[10]),np.array(bigenergy[10]),label=labels[10], lw=0.5)
#ax.legend(handles=lines,fontsize=8, loc=2)#,line2,line3,line4,line5,line6, line7, line8, line9, line10], fontsize=8, loc=2)
ax.legend(labels, fontsize=12, loc=2)
ax.set_xlabel(r"time, $t_{\mbox{mpcd}}$",fontsize=FS)
ax.set_ylabel(r"Total energy",fontsize=FS)
ax.set_title(r"System size 100x100, varying thermostating techniques",fontsize=FS)
plt.savefig(cwd+"/100x100thermostat.pdf",format='pdf', dpi = 'figure')
plt.savefig(cwd+"/100x100thermostat.png",format='png', dpi = 'figure')
plt.close(fig)