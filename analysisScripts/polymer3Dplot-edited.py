import matplotlib
# matplotlib.use('Agg')
from pylab import *
from numpy import ma
from subprocess import call
from mpl_toolkits.mplot3d import axes3d
import os
import sys

###########################################################
### Plots 3D polymer conformation and director field
###########################################################

###########################################################
### Read arguments
###########################################################
xyzSize=zeros( 3,dtype=int )
print( "Arguments:" )
for arg in sys.argv:
	print( "\t" + arg )
mpcdDataPath = sys.argv[1]	# Path to directory containing MPCD data
###########################################################
### Initialize
###########################################################
c = 0.5		#Length of director lines approx 0.5
keepFrames=0	#0=don't keep (delete) frames; 1=keep frames
start = 0		# Show after this number (in MPCD)
finish = 999999		# Finish after this number
qx = 1.0		# Only show every qx arrow in x
qy = 1.0	# Only show every qy arrow in y
qz = 1.0		# Only show every qz arrow in z

###########################################################
FS =25
TLS = 20		# Tick label size
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
#############################################################
#Animation stuff
bitrate=5000
framerate=4		#Number of frames per second in the output video
codec='libx264'		#Other options include mpeg4
suffix='.mp4'
###########################################################
### Read input files
###########################################################
# print( '\tReading md input md.inp ...' )
print( '\tReading md input md.inp ...' )
file=mpcdDataPath+'/md.inp'
if not os.path.isfile(file):
	print("%s not found."%file)
	exit()
file=open(file,"r")
for i in range(11):
	buff=file.readline()
buff=file.readline().split('=')[-1].replace('(','').replace(')','').split(',')
for i in range(3):
	xyzSize[i] = int(buff[i])
	if (xyzSize[i] == 0) :
		xyzSize[i] += 1
for i in range(5):
	buff=file.readline()
buff=buff.split('=')
dtMD=float(buff[-1])
for i in range(59):
	buff=file.readline()
buff=file.readline().split()
monoN=int(buff[-1].replace('(','').replace(')',''))
monoNF=float(monoN)
for i in range(36):
	buff=file.readline()
buff=file.readline().split(',')[-1].split(')')
outMD=int(buff[0])
file.close()

MDMPCD=int(0.1/dtMD)                             # number of MD timesteps per MPCD timestep
dirSOut=int(outMD/MDMPCD) 

nemName=mpcdDataPath+'/directorfield.dat'
# print( '\tFinding vmd.vtf ...' )
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
# Data
XYZ = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
DIR = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
S = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
POLY = zeros(shape=(3,monoN),dtype=float)

# Figure
fig1 = plt.figure(1)
ax = fig1.add_subplot(111, projection='3d')
ax = fig1.gca(projection='3d')

#Create the colorbar
CB = ax.scatter3D([0,0], [0,0], zs=[0,0], c=[0,1], s=0, cmap=myMap)
cb = fig1.colorbar(CB)
cb.ax.set_ylabel(r'$S$', fontsize = FS)

###########################################################
### Read the data for animation
###########################################################
labX='x'
labY='y'
labZ='z'

# print( 'Read file for figures ...' )
nemFile = open(nemName,"r")
# print( '\tToss headers ...' )
for i in range(13):
	#toss header
	line = nemFile.readline()
	#print line
polyFile = open(polyName,"r")
while True:
	line = polyFile.readline()
	if(line[0]=='#'):
		pass
	else:
		break
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
			XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx) + 0.5
			XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy) + 0.5
			XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz) + 0.5
			DIR[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
			DIR[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
			DIR[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)
			S[int(Qx)][int(Qy)][int(Qz)] = float(s)

	if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
		# Read polymer
		line = polyFile.readline()
		line = polyFile.readline()
		for k in range(monoN):
			line = polyFile.readline()
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
		j=j+1
		if j>finish:
			break

		#Save the instantaneous or current velocity field frame
		# Save frame
		n=n+1
		print( "\t\tPlot frame %d"%(n) )
		fig1 = plt.figure(1)
		plt.cla()
		for x in range(xyzSize[0]):
			for y in range(xyzSize[1]):
				for z in range(xyzSize[2]):
					if( x%qx==0 and y%qy==0 and z%qz==0 ):
						if( S[x][y][z]>0.05 ):
							xline=[ XYZ[0][x][y][z]-c*DIR[0][x][y][z],XYZ[0][x][y][z]+c*DIR[0][x][y][z] ]
							yline=[ XYZ[1][x][y][z]-c*DIR[1][x][y][z],XYZ[1][x][y][z]+c*DIR[1][x][y][z] ]
							zline=[ XYZ[2][x][y][z]-c*DIR[2][x][y][z],XYZ[2][x][y][z]+c*DIR[2][x][y][z] ]
							plot( xline,yline,zline,color=[240./255, 191./255., 79./255.],alpha=0.5 )
		for k in range(monoN-1):
			if(fabs(POLY[0][k]-POLY[0][k+1])<0.5*xyzSize[0] and fabs(POLY[1][k]-POLY[1][k+1])<0.5*xyzSize[1] and fabs(POLY[2][k]-POLY[2][k+1])<0.5*xyzSize[2]):
				plot( [POLY[0][k],POLY[0][k+1]],[POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='k',linewidth=3 )
		ax.set_xlabel(r'$%s$'%labX, fontsize = FS)
		ax.set_ylabel(r'$%s$'%labY, fontsize = FS)
		ax.set_zlabel(r'$%s$'%labZ, fontsize = FS)
		ax.set_xlim3d(0,xyzSize[0])
		ax.set_ylim3d(0,xyzSize[1])
		ax.set_zlim3d(0,xyzSize[2])
		ax.grid(False)
		plt.axis(xmax=xyzSize[0], xmin=0, ymax=xyzSize[1], ymin=0) #plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
		name='frame%04d.png'%(n)
		savefig( name, format='png', dpi=600 )

		#Zero matrix
		DIR= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
		S = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
		i=0
		# print( "\tRead frame %d"%(n) )
nemFile.close()
polyFile.close()

###########################################################
### Animate data
###########################################################
# print( "\tAnimating ..." )
name='3DorientationPolymer_animation%s'%suffix
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
	myCommand="rm frame*.png"
	call(myCommand,shell=True)
