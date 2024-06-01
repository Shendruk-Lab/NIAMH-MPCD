from matplotlib import pyplot as plt
from pylab import *
from numpy import ma
from subprocess import call
import os
import sys
rcParams["figure.autolayout"] = True
###########################################################
### Plots 2D averaging over user defined direction
###########################################################

###########################################################
### Read arguments
###########################################################

xyzSize=zeros( 3,dtype=int )
print( "Arguments:" )
for arg in sys.argv:
	print( "\t" + arg )
mpcdDataPath = sys.argv[1]		# Path to directory containing MPCD data (input.json, md.inp, directorfield.dat, folder for vmd file)
###########################################################
### Initialize
###########################################################
projection = 'z'		# The direction to project over
c = 0.5		#Length of director lines approx 0.5
keepFrames=1	#0=don't keep (delete) frames; 1=keep frames
start = 0		# Show after this number (in MPCD)
finish = 999999999		# Finish after this number

############################################################
qx = 0.5
qy = 0.5

###########################################################
### Format and style
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
myMap=ed.plasma
# Adjust line width
myLW=1.0


FS =25
TLS = 20		# Tick label size
#Animation stuff
bitrate=5000
framerate=4		#Number of frames per second in the output video
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

###########################################################
### Read input files
###########################################################
print( '\tReading md input md.inp ...' )
file=mpcdDataPath+'/md.inp'
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
for i in range(5):
	buff=file.readline()
	# reads params lattice to dt
buff=buff.split('=')
dtMD=float(buff[-1])
for i in range(59):
	buff=file.readline()
	# reads remaining lines
buff=file.readline().split()
monoN=int(buff[-1].replace('(','').replace(')',''))
monoNF=float(monoN)
for i in range(36):
	buff=file.readline()
buff=file.readline().split(',')[-1].split(')')
outMD=int(buff[0])
file.close()

MDMPCD=int(0.1/dtMD)                             # number of MD timesteps per MPCD timestep
dirSOut=int(outMD/MDMPCD)        				 # the frequency of mpcd output 

nemName=mpcdDataPath+'/directorfield.dat'
print( '\tFinding vmd.vtf ...' )
check=0
t=0
for r, d, f in os.walk(mpcdDataPath):
	for file in f:
		if '-vmd.vtf' in file:
			check+=1
			polyName = os.path.join(r,file)
if(not check):
	print("VMD file not found.")
	exit()

# Number of MD time steps outputted per MPCD outputs
outMD_in_MPCD_iterations=outMD/MDMPCD
MDperMPCD=dirSOut/outMD_in_MPCD_iterations
###########################################################
### Initialize
###########################################################
XYZ = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
DIR = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
# S = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)

## for plotting scalar field using imshow
S = zeros(shape=(xyzSize[0],xyzSize[1]),dtype=float)

if(projection=='x' or projection=='X'):
	projDim=0
	labX='y'
	labY='z'
	xySize=array([xyzSize[1],xyzSize[2]])
elif(projection=='y' or projection=='Y'):
	projDim=1
	labX='x'
	labY='z'
	xySize=array([xyzSize[0],xyzSize[2]])
elif(projection=='z' or projection=='Z'):
	projDim=2
	labX='x'
	labY='y'
	xySize=array([xyzSize[0],xyzSize[1]])
else:
	print("Error: Projection direction must be x, y or z")
	exit()
POLY = zeros(shape=(3,monoN),dtype=float)

# Setup figure object
tight_layout()
fig,ax = plt.subplots()
#Create the colorbar
CS3 = imshow(zeros(shape=(xySize[0],xySize[1]),dtype=float), cmap=myMap, vmin=0, vmax=1, extent=[0,xySize[0],0,xySize[1]])
cb = fig.colorbar(CS3)
cb.ax.set_ylabel(r'$S$', fontsize = FS)
###########################################################
### Read the data for animation
###########################################################
# print( 'Read file for figures ...' )
nemFile = open(nemName,"r")
# print( '\tToss headers ...' )
for i in range(13):
	#toss header
	line = nemFile.readline()
	#print line
polyFile = open(polyName,"r")
j=0
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

# print( '\tRead data ...' )
i=0
j=0
n=-1
# print( "\tRead frame %d"%(n+1) )
while nemFile:
	i=i+1
	# Read nematic field
	line = nemFile.readline()
	if( not line ):
		break
	else:
		if j>=start:
			t,Qx,Qy,Qz,Vx,Vy,Vz,s = line.split("\t",8)
			XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx)
			XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy)
			XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz)
			DIR[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
			DIR[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
			DIR[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)
			# S[int(Qx)][int(Qy)][int(Qz)] = float(s)
			S[int(Qx)][int(Qy)] = float(s)

	if i==xyzSize[0]*xyzSize[1]*xyzSize[2]: #because polymer only read once per sweep of field
		# Read polymer
		line = polyFile.readline()
		line = polyFile.readline()
		# print(line)
		for k in range(monoN):
			line = polyFile.readline()
			# print line
			t,Qx,Qy,Qz = line.split()
			POLY[0][k]=float(Qx)
			POLY[1][k]=float(Qy)
			POLY[2][k]=float(Qz)
		# Throw away the next MDperMPCD-1 frames
		for l in range(int(MDperMPCD-1)):
			line = polyFile.readline()
			line = polyFile.readline()
			for k in range(monoN):
				line = polyFile.readline()
		# Calculate the cm
		cm=zeros(shape=(3),dtype=float)
		for d in range(3):
			for i in range(monoN):
				cm[d]+=POLY[d][i]
			cm[d]/=float(monoN)
		j=j+1
		if j>finish:
			break

		# Save the instantaneous or current velocity field frame
		# Save frame
		n=n+1
		print( "\t\tPlot frame %d"%(n) )
		plt.cla()
		if(projDim==0):
			x=int(cm[0])
			for y in range(xyzSize[1]):
				for z in range(xyzSize[2]):
					if( S[x][y][z]>0.05 ):
						yline=[ XYZ[1][x][y][z]-c*DIR[1][x][y][z],XYZ[1][x][y][z]+c*DIR[1][x][y][z] ]
						zline=[ XYZ[2][x][y][z]-c*DIR[2][x][y][z],XYZ[2][x][y][z]+c*DIR[2][x][y][z] ]
						plot( yline,zline,color=myMap(S[x][y][z],1),linewidth=myLW )
			for k in range(monoN-1):
				if(fabs(POLY[1][k]-POLY[1][k+1])<0.5*xyzSize[1] and fabs(POLY[2][k]-POLY[2][k+1])<0.5*xyzSize[2]):
					plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='k',linewidth=2 )
		elif(projDim==1):
			y=int(cm[1])
			for x in range(xyzSize[0]):
				for z in range(xyzSize[2]):
					if( S[x][y][z]>0.05 ):
						xline=[ XYZ[0][x][y][z]-c*DIR[0][x][y][z],XYZ[0][x][y][z]+c*DIR[0][x][y][z] ]
						zline=[ XYZ[2][x][y][z]-c*DIR[2][x][y][z],XYZ[2][x][y][z]+c*DIR[2][x][y][z] ]
						plot( xline,zline,color=myMap(S[x][y][z],1),linewidth=myLW )
			for k in range(monoN-1):
				if(fabs(POLY[0][k]-POLY[0][k+1])<0.5*xyzSize[0] and fabs(POLY[2][k]-POLY[2][k+1])<0.5*xyzSize[2]):
					plot( [POLY[0][k],POLY[0][k+1]],[POLY[2][k],POLY[2][k+1]],color='k',linewidth=2 )
		elif(projDim==2):
			z=int(cm[2])
			for x in range(xyzSize[0]):
				for y in range(xyzSize[1]):
					# if( S[x][y][z]>0.05 ):
					if( S[x][y]>0.05 ):
						xline=[ XYZ[0][x][y][z]-c*DIR[0][x][y][z],XYZ[0][x][y][z]+c*DIR[0][x][y][z] ]
						yline=[ XYZ[1][x][y][z]-c*DIR[1][x][y][z],XYZ[1][x][y][z]+c*DIR[1][x][y][z] ]
						# plot( xline,yline,color=myMap(S[x][y][z],1),linewidth=myLW )
						plot( xline,yline,color="white",linewidth=myLW )
			imshow(S.T,cmap=myMap,origin='lower',vmin=0,vmax=1)
			# circle = Circle((cm[j-1][0],cm[j-1][1]),radius = 0.5, facecolor = "none", edgecolor = "black",zorder =1)
			# if (j-1>0):
			# 	for k in range(j-1):
			# 		plot( [cm[k][0],cm[k+1][0]],[cm[k][1],cm[k+1][1]],color='green',linewidth=2,zorder=2 )
			for k in range(monoN-1):
				if(fabs(POLY[0][k]-POLY[0][k+1])<0.5*xyzSize[0] and fabs(POLY[1][k]-POLY[1][k+1])<0.5*xyzSize[1]):
					plot( [POLY[0][k],POLY[0][k+1]],[POLY[1][k],POLY[1][k+1]],color='black',linewidth=2,zorder=2 )


		# axis(xmax=xySize[0], xmin=0, ymax=xySize[1], ymin=0)
		# plt.xticks(ma.arange(0, xySize[0]+1, 5))
		# plt.yticks(ma.arange(0, xySize[1]+1, 5))
		# ax.set_yticklabels([])
		# ax.set_xticklabels([])
		# ax.set_xlabel(r'$%s$'%labX, fontsize = FS)
		# ax.set_ylabel(r'$%s$'%labY, fontsize = FS)
		# ax.tick_params(axis=u'both', which=u'both',length=0)
		# ax.set_aspect('equal', adjustable='box')


		# ax = plt.gca()
		# ax.add_patch(circle)

		axis('off')
		grid(False)
		subplots_adjust(left=0, bottom=0, right=1, top=1)
		name='frame%04d.png'%(n)
		savefig( name, format='png', dpi=600, bbox_inches='tight')
		i=0

nemFile.close()
polyFile.close()

###########################################################
### Animate data
###########################################################
# print( "\tAnimating ..." )
name='2DPolymerProj_animation%s'%suffix
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
	myCommand="rm frame*.png"
	call(myCommand,shell=True)
