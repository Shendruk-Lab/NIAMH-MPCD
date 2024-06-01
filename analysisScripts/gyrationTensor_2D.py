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
### Calculates radius of gyration measures for 2D
###########################################################

###########################################################
### Read arguments
###########################################################
def read(cwd):
	mpcdDataPath = cwd
	file=mpcdDataPath+'/md.inp'
	if not os.path.isfile(file):
		print("%s not found."%file)
		exit()
	file=open(file,"r")
	#for i in range(16):
	#	buff=file.readline()
	#buff=file.readline().split('=')
	#dtMD=float(buff[-1])
	for i in range(29): #was 21/38
		buff = file.readline()
	buff = file.readline().split()
	print(buff)
	bend = float(buff[-1].replace('(','').replace(')',''))
	for i in range(8): #was 21/38
		buff = file.readline()
	buff = file.readline().split()
	act = float(buff[-1].replace('(','').replace(')',''))
	buff = file.readline().split()
	chunks = int(buff[-1].replace('(','').replace(')',''))
	file.close()
 
	return bend, act, chunks

def func(cwd):
	xyz=np.zeros( 2,dtype=int )
	#print( "Arguments:" )
	#for arg in sys.argv:
	#	print( "\t" + arg )
	#mpcdDataPath = sys.argv[1]	# Path to directory containing md.inp
	#mdDataPath = sys.argv[2]	# Path to directory containing vmd.vtf file
	#xyz[0] = int(sys.argv[3])	# System size
	#xyz[1] = int(sys.argv[4])	# System size
	mpcdDataPath = cwd
	mdDataPath = os.getcwd()
	print(mdDataPath)
	xyz[0] = 100
	xyz[1] = 100

	###########################################################
	### Initialize
	###########################################################
	_x=0
	_y=1
	###########################################################
	### Read the data
	###########################################################
	# print( '\tReading md input md.inp ...' )
	file=mpcdDataPath+'/md.inp'
	if not os.path.isfile(file):
		print("%s not found."%file)
		exit()
	file=open(file,"r")
	for i in range(16):
		buff=file.readline()
	buff=file.readline().split('=')
	dtMD=float(buff[-1])
	for i in range(59):
		buff=file.readline()
	buff=file.readline().split()
	numMono=int(buff[-1].replace('(','').replace(')',''))
	print(numMono)
	numMonoF=float(numMono)
	for i in range(36):
		buff=file.readline()
	buff=file.readline().split(',')[-1].split(')')
	outMD=int(buff[0])
	file.close()
	dtMD_out=dtMD*outMD

	# print( '\tReading vmd.vtf ...' )
	timeMD=[]
	posWrap=[]
	check=0
	t=0
	for r, d, f in os.walk(mdDataPath):
		for file in f:
			if '-vmd.vtf' in file:
				check+=1
				print(file)
				file=os.path.join(r,file)
				file=open(file,"r")
				for i in range(77): #was 75 but not working
					buff=file.readline()
				while file:
					buff=file.readline()
					if( not buff ):
						break
					timeMD.append(t)	#In MPCD time units
					posWrap.append(np.zeros(shape=(numMono,2),dtype=float))
					for i in range(numMono):
						buff=file.readline().split()
						#print(buff)
						for d in range(2):
						#print(t, i, d)
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

	# print( "\tCalculating CM ..." )
	cm=np.zeros(shape=(2,mdSteps),dtype=float)
	for t in range(mdSteps):
		for d in range(2):
			for i in range(numMono):
				cm[d][t]+=pos[i][d][t]
			cm[d][t]/=numMonoF
	
	# print( "\tCalculating gyration tensor ..." )
	rg=np.zeros(mdSteps,dtype=float)
	rg2Tensor=np.zeros(shape=(mdSteps,2,2),dtype=float)
	rg_0=np.zeros(mdSteps,dtype=float)
	rg_1=np.zeros(mdSteps,dtype=float)
	RelShapeAnisotropy=np.zeros(mdSteps,dtype=float)
	for t in range(mdSteps):
		rg2=0.0
		for i in range(numMono):
			rg2T=np.zeros(shape=(2,2),dtype=float)
			for j in range(2):
				for k in range(2):
					rg2T[j][k]=(pos[i][j][t]-cm[j][t])*(pos[i][k][t]-cm[k][t])
					rg2Tensor[t][j][k]+=rg2T[j][k]
			rg2+=rg2T[_x][_x]+rg2T[_y][_y]
		rg[t]=np.sqrt( rg2/numMonoF )
		for j in range(2):
			for k in range(2):
				rg2Tensor[t][j][k]/=numMonoF
		w2, v = LA.eig(rg2Tensor[t])
		#Sometimes eig() gets the order wrong so just make sure
		w2[::-1].sort()	
		if w2[0]>=0.0:
			rg_0[t]=np.sqrt(w2[0])
		if w2[1]>=0.0:
			rg_1[t]=np.sqrt(w2[1])
		RelShapeAnisotropy[t]= math.pow((w2[0]-w2[1]),2)/math.pow((w2[0]+w2[1]),2)

	# print( "\tCalculating end-to-end vector ..." )
	ree=np.zeros(mdSteps,dtype=float)
	ree_v=np.zeros(shape=(mdSteps,2),dtype=float)
	for t in range(mdSteps):
		for i in range(2):
			ree_v[t][i]=pos[numMono-1][i][t]-pos[0][i][t]
			ree[t]+=ree_v[t][i]*ree_v[t][i]
		ree[t]=np.sqrt(ree[t])

	##creating a file which includes results
	outfile = open(mdDataPath+"/gyration.dat", 'w')
	outfile.write("Gyration tensor measures and end to end distance with time \n")
	outfile.write("time\trg\trg_0\trg_1\tk^2\tree\n")
	for t in range(mdSteps):
		outfile.write("%e\t%e\t%e\t%e\t%e\t%e\n"%(timeMD[t],rg[t], rg_0[t], rg_1[t],RelShapeAnisotropy[t],ree[t]))	
	# close output file
	outfile.close()

def cmfunc(cwd,prevwd):
	xyz=np.zeros( 2,dtype=int )
	#print( "Arguments:" )
	#for arg in sys.argv:
	#	print( "\t" + arg )
	#mpcdDataPath = sys.argv[1]	# Path to directory containing md.inp
	#mdDataPath = sys.argv[2]	# Path to directory containing vmd.vtf file
	#xyz[0] = int(sys.argv[3])	# System size
	#xyz[1] = int(sys.argv[4])	# System size
	mpcdDataPath = prevwd
	mdDataPath = cwd #os.getcwd()
	print(mdDataPath)
	xyz[0] = 100
	xyz[1] = 100

	###########################################################
	### Initialize
	###########################################################
	_x=0
	_y=1
	###########################################################
	### Read the data
	###########################################################
	# print( '\tReading md input md.inp ...' )
	file=mpcdDataPath+'/md.inp'
	if not os.path.isfile(file):
		print("%s not found."%file)
		exit()
	file=open(file,"r")
	for i in range(16):
		buff=file.readline()
	buff=file.readline().split('=')
	dtMD=float(buff[-1])
	for i in range(59):
		buff=file.readline()
	buff=file.readline().split()
	numMono=int(buff[-1].replace('(','').replace(')',''))
	print(numMono)
	numMonoF=float(numMono)
	for i in range(36):
		buff=file.readline()
	buff=file.readline().split(',')[-1].split(')')
	outMD=int(buff[0])
	file.close()
	dtMD_out=dtMD*outMD
	print(dtMD_out, " this should be 25")

	# print( '\tReading vmd.vtf ...' )
	timeMD=[]
	posWrap=[]
	check=0
	t=0
	for r, d, f in os.walk(mdDataPath):
		for file in f:
			if '-vmd.vtf' in file:
				check+=1
				print(file)
				file=os.path.join(r,file)
				file=open(file,"r")
				for i in range(77): #was 75 but not working
					buff=file.readline()
				while file:
					buff=file.readline()
					if( not buff ):
						break
					timeMD.append(t)	#In MPCD time units
					posWrap.append(np.zeros(shape=(numMono,2),dtype=float))
					for i in range(numMono):
						buff=file.readline().split()
						#print(buff)
						for d in range(2):
						#print(t, i, d)
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

	# print( "\tCalculating CM ..." )
	cm=np.zeros(shape=(2,mdSteps),dtype=float)
	print(np.shape(cm))
	for t in range(mdSteps):
		for d in range(2):
			for i in range(numMono):
				cm[d][t]+=pos[i][d][t]
			cm[d][t]/=numMonoF
	
	return cm, dtMD_out

#func(os.getcwd())
