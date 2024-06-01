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
def func(mdDataPath, xyz0, xyz1, numMono):
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


	# Calculate stats
	# averaging over backbone and a function of time
	avK = zeros(shape=(mdSteps),dtype=float)
	stdK = zeros(shape=(mdSteps),dtype=float)
	cnts = zeros(shape=(mdSteps),dtype=int)
	for t in range(mdSteps):
		for i in range(1,numMono-1):
			cnts[t]+=1.
			avK[t]+=curvature[t][i]
			stdK[t]+=curvature[t][i]*curvature[t][i]
		if cnts[t]>2:
			stdK[t]=sqrt((cnts[t]*stdK[t]-avK[t]*avK[t])/(cnts[t]*(cnts[t]-1.)))
		else:
			stdK[t]=0.0
		if cnts[t]>1:
			avK[t]/=cnts[t]
	##creating a file which includes results
	outfile = open(mdDataPath+"/curvature-time.dat", 'w')
	outfile.write("Average curvature with time \n")
	outfile.write("time\tcurvature\tstd\n")
	for t in range(mdSteps):
		outfile.write("%e\t%e\t%e\n"%(t,avK[t],stdK[t]))	
	outfile.close()

	time = [i*20 for i in range(mdSteps)]
	fig,ax = plt.subplots(1)
	ed.errorbar_fill(np.array(time),np.array(avK),yerr=np.array(stdK))
	# ax.plot(time,array)
	ax.set_xlabel(r"time, $t_{\mbox{mpcd}}$",fontsize=FS)
	ax.set_ylabel(r"Curvature, $\kappa$",fontsize=FS)
	plt.savefig(mdDataPath+"/curvature-time.pdf",format='pdf', dpi = 'figure')
	plt.close(fig)
	#############################################################################
	# averaging over time and a function of position along backbone
	#data = np.genfromtxt(warmupPath+'warmup1.dat', skip_header=1)
	#warmup = data.T
	#warmup = int(warmup)
	warmup=0
	avK = [mean(x) for x in zip(*curvature[warmup:])]
	stdK = [std(x) for x in zip(*curvature[warmup:])]
	##creating a file which includes results
	outfile = open(mdDataPath+"/curvature-L.dat", 'w')
	outfile.write("Average curvature along backbone \n")
	outfile.write("position\tcurvature\tstd\n")
	for i in range(numMono):
		outfile.write("%e\t%e\t%e\n"%(i,avK[i],stdK[i]))	
	outfile.close()

	L = [i for i in range(numMono)]
	fig,ax = plt.subplots(1)
	ed.errorbar_fill(np.array(L),np.array(avK),yerr=np.array(stdK))
	# ax.plot(L,array)
	ax.set_xlabel(r"Monomer index, $m$",fontsize=FS)
	ax.set_ylabel(r"Curvature, $\kappa$",fontsize=FS)
	plt.savefig(mdDataPath+"/curvature-L.pdf",format='pdf', dpi = 'figure')
	plt.close(fig)