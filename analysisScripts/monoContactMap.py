from matplotlib import pyplot as plt
import matplotlib.colors as colors
from pylab import *
from numpy import ma
from subprocess import call
import os
import sys
#import skimage
rcParams["figure.autolayout"] = True
FS =15
###########################################################
### Plots 2D averaging over user defined direction
###########################################################

###########################################################
### Read arguments
###########################################################
def func(mdinpDataPath, bend, act, chunks, log):
	xyzSize=zeros( 3,dtype=int )
	print( "Arguments:" )
	for arg in sys.argv:
		print( "\t" + arg )
	#mdinpDataPath = sys.argv[1]		# Path to directory containing mdinp data (input.json, md.inp)
	#act = sys.argv[2]
	#chunks = sys.argv[3]
	if log==True:
		thresh = 2
	else:
		thresh = 3		# HAVE BEEN CHANGING TO 2 FOR LOGARITHMIC TO MAKE LESS NOISY CAREFUL!!!

	###########################################################
	### Format and style
	###########################################################
	# Use our custom style
	plt.style.use('shendrukGroupStyle')
	# Use our custom colours
	import shendrukGroupFormat as ed
	myMap=ed.viridis

	###########################################################
	### Read input files
	###########################################################
	print( '\tReading md input md.inp ...' )
	file=mdinpDataPath+'/md.inp'
	if not os.path.isfile(file):
		print("%s not found."%file)
		exit()
	file=open(file,"r")
	for i in range(11):
		buff=file.readline()
		# reads all lines up to after randomSeed
	buff=file.readline().split('=')[-1].replace('(','').replace(')','').split(',')
	# reads dimensions of system
	for i in range(3):
		xyzSize[i] = int(buff[i])
		if (xyzSize[i] == 0) :
			xyzSize[i] += 1
	#for i in range(5):
	#	buff=file.readline()
	#	# reads params lattice to dt
	#buff=buff.split('=')
	#dtMD=float(buff[-1])	#DONT NEED
	#for i in range(21): 
	#    buff = file.readline()
	#buff = file.readline().split() #for the chunk colours
	#dStrength=float(buff[-1].replace('(','').replace(')',''))	#DONT NEED
	#buff = file.readline().split()
	#dChunks = int(buff[-1].replace('(','').replace(')',''))
	for i in range(64): #36+23+6
		buff=file.readline()
		print(buff)
		# reads remaining lines
	buff=file.readline().split()
	monoN=int(buff[-1].replace('(','').replace(')','')) #NEED !!
	monoNF=float(monoN)	#DONT NEED
	for i in range(36):
		buff=file.readline()
	buff=file.readline().split(',')[-1].split(')')
	outMD=int(buff[0]) #PROBS DONT NEED
	file.close()

	def Threshold(p1,p2):
		"""
		For calculating distance between monomers and determining if they are
		in 'contact' according to the threshold set.
		Args:
			p1 (_type_): position array of one particle
			p2 (_type_): same but for second particle
		Returns:
			1 for yes or 0 for no
		"""
		d = norm(p2-p1)
		if d<thresh:
			return 1
		else:
			return 0
		
	map = zeros(shape=(monoN,monoN),dtype=int)
	for i in range(20):
		print(i)
		repeat=i+1
		polydatapath = str(mdinpDataPath+"/result"+str(repeat))
		print(polydatapath)

		print( '\tFinding vmd.vtf ...' )
		check=0
		t=0
		for r, d, f in os.walk(polydatapath):
			for file in f:
				if '-vmd.vtf' in file:
					check+=1
					polyName = os.path.join(r,file)
		if(not check):
			print("VMD file not found.")
			exit()
		POLY = zeros(shape=(monoN,3),dtype=float)
		polyFile = open(polyName,"r")
		#j=0
		while True:
			line = polyFile.readline()
			# print(j)
			if(line[0]=='#'):
				pass
			else:
				break
			# j=j+1
		line = polyFile.readline()
		line = polyFile.readline()
		#print(line)
		j=-1
		while j<200:
			j+=1
			# zero temporary map
			tempmap = zeros(shape=(monoN,monoN),dtype=int)
			# Read polymer
			line = polyFile.readline()
			line = polyFile.readline()
			#print(j)
			for k in range(monoN):
				line = polyFile.readline()
				# print line
				#print(line)
				t,Qx,Qy,Qz = line.split()
				POLY[k][0]=float(Qx)
				POLY[k][1]=float(Qy)
				POLY[k][2]=float(Qz)
			for x in range(monoN):
				for y in range(monoN):
					p1 = POLY[x]
					p2 = POLY[y]
					if not x==y:
						if not x==(y-1)%monoN:
							if not x==(y+1)%monoN:
								tempmap[x][y]=Threshold(p1,p2)
			map+=tempmap
		polyFile.close()

	prob=map/norm(map)
	# Setup figure object
	tight_layout()
	fig,ax = plt.subplots()
	plt.cla()
	if log==False:
		CS3 = imshow(prob.T, cmap=myMap,origin='lower')
		cb = fig.colorbar(CS3)
		cb.ax.set_ylabel(r'Probablitiy of monomers being in contact', fontsize = FS)
	else:
		CS3 = ax.pcolor(prob.T, norm=colors.LogNorm(), cmap=myMap)
		cb.ax.set_ylabel(r'Probablitiy of monomers being in contact (logarithmic)', fontsize = FS)
		cb.ax.set_ylabel(r'Probablitiy of monomers being in contact', fontsize = FS)
	axis('off')
	grid(False)
	subplots_adjust(left=0, bottom=0, right=1, top=1)
	name=str('/home/s1954660/Desktop/summer/analysis/LOGcontactmaps/LOGcontactmap-'+str(act)+'-'+str(chunks)+'.png')
	savefig( name, format='png', dpi=600, bbox_inches='tight')