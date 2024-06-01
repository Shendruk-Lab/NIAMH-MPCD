import matplotlib
from pylab import *
from numpy import ma
from subprocess import call
from numpy import linalg as LA
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors
import os
import csv
import sys
import numpy as np
import math
###########################################################
### Calculates and plot curvature for 2D
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
myMap=ed.plasma
# Adjust line width
myLW=1.0
FS=20
###########################################################
### Read arguments
###########################################################
def func(mdDataPath, xyz0, xyz1, numMono,log):
    xyz=np.zeros( 2,dtype=int )
	#print( "Arguments:" )
	#for arg in sys.argv:
	#	print( "\t" + arg )
	#mdDataPath = sys.argv[1]	# Path to directory containing vmd.vtf file
	#warmupPath = sys.argv[2]	# Path to directory containing warmup time file
	#xyz[0] = int(sys.argv[3])	# System size
	#xyz[1] = int(sys.argv[4])	# System size
    xyz[0] = xyz0
    xyz[1] = xyz1
	#numMono = int(sys.argv[5])
	###########################################################
	### Read the data
	###########################################################
	# print( '\tReading vmd.vtf ...' )
    timeMD=[]
    posWrap=[]
    check=0
    t=0
    for r, d, f in os.walk(mdDataPath):
        for file in f:
            if '-vmd.vtf' in file:
                check+=1
                file=os.path.join(r,file)
                file=open(file,"r")
                for i in range(77): #was 75
                    buff=file.readline()
                while file:
                    buff=file.readline()
                    if( not buff ):
                        break
                    timeMD.append(t)	#In MPCD time units
                    posWrap.append(np.zeros(shape=(numMono,2),dtype=float))
                    for i in range(numMono):
                        buff=file.readline().split()
                        for d in range(2):
                            posWrap[t][i][d]=float(buff[1+d])
                    t+=1
                    buff=file.readline()
                file.close()
    if(not check):
        print("VMD file not found.")
        exit()

	###########################################################
	### Analysis
	###########################################################
    mdSteps=len(timeMD)

	# print( "\tUnwrap positions ..." )
    pos=np.zeros(shape=(numMono,2,mdSteps),dtype=float)
    for t in range(mdSteps):
        for d in range(2):
            pos[0][d][t]=posWrap[t][0][d]
        for i in range(1,numMono):
            for d in range(2):
                dx=posWrap[t][i][d]-posWrap[t][i-1][d]
                if(dx>0.5*xyz[d]):
                    dx-=xyz[d]
                elif(dx<-0.5*xyz[d]):
                    dx+=xyz[d]
                pos[i][d][t]=pos[i-1][d][t]+dx


    sep = zeros(shape=(numMono-1),dtype=float)
    unitTang = zeros(shape=(2,numMono-1),dtype=float)
    curvature = zeros(shape=(mdSteps,numMono),dtype=float)
	# Calculate the unit tangent vectors from forward differences (from 0 to numMono-1, i.e. none for the last monomer)
    for t in range(mdSteps):
        for i in range(numMono-1):
            sep[i]=0.0
            for d in range(2):
                dx=pos[i+1][d][t]-pos[i][d][t]
                unitTang[d][i]=dx
                sep[i]+=dx*dx
            sep[i]=sqrt(sep[i])
            for d in range(2):
                unitTang[d][i]/=sep[i]
		# Calculate the curvature (assume the ends have zero curvature)
        for i in range(1,numMono-1):
            dtds=zeros(shape=(2),dtype=float)
            for d in range(2):
                dtds[d]=(unitTang[d][i]-unitTang[d][i-1])/sep[i-1]
            curvature[t][i]=sqrt(dtds[0]*dtds[0]+dtds[1]*dtds[1])

	##creating a file which includes results
    outfile = open(mdDataPath+"/chemo-curv.dat", 'w')
    outfile.write("Chemograph of curvature, time and index. \n")
    outfile.write("time\tindex\tcurvature\n")
    for t in range(mdSteps):
        for i in range(numMono):
            outfile.write("%e\t%e\t%e\n"%(t,i,curvature[t][i]))	
    outfile.close()

    fig,ax = plt.subplots(1)
    plt.cla()
    if log==True:
        chemograph = ax.pcolor(curvature,norm=colors.LogNorm(),cmap=myMap)
    else:
        chemograph = ax.pcolor(curvature,cmap=myMap)
    #chemograph = imshow(curvature,cmap=myMap, origin='lower', aspect='auto')
    cb=fig.colorbar(chemograph)
    #cb.ax.set_ylabel(r"Curvature, $\kappa$",fontsize=FS)
    if log==True:
        cb.ax.set_ylabel(r"Curvature, $\kappa$ - LOGARITHMIC",fontsize=FS)
    else:
        cb.ax.set_ylabel(r"Curvature, $\kappa$",fontsize=FS)
    xlabel(r"Monomer index",fontsize=FS)
    ylabel(r"time, $t_{\mbox{mpcd}}$/500",fontsize=FS)
    if log==True:
        plt.savefig(mdDataPath+"/LOGchemo_curv.pdf",format='pdf', dpi = 'figure')
    else:
        plt.savefig(mdDataPath+"/chemo_curv.pdf",format='pdf', dpi = 'figure')
    plt.close(fig)
 
#log=False   
#func(os.getcwd(), 100, 100, 24,log)
