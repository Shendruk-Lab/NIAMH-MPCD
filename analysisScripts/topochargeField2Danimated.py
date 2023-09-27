from pylab import *
from numpy import ma
from subprocess import call
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import json

###########################################################
### Plots 2D averaging over user defined direction
###########################################################

###########################################################
### Read arguments
###########################################################
print( "Arguments:" )
for arg in sys.argv:
	print( "\t" + arg )
if(len(sys.argv) != 13): #check for correct number of arguments to be idiot proof
	print("Error: This script expects 12 arguments.\nTerminating\n")
	quit()
directorData = sys.argv[1]		# Name of the field data
topoData = sys.argv[2]			# Name of the topo data
inputName = sys.argv[3]			# Input json file to read inputs
start = int(sys.argv[4])		# Average after this number
finish = int(sys.argv[5])		# Average before this number
qx = int(sys.argv[6])		# Only show every qx arrow in x
qy = int(sys.argv[7])		# Only show every qy arrow in y
avdim = sys.argv[8]			# Dimension to average over
c = float(sys.argv[9])		#Length of director lines approx 0.5
c1 = float(sys.argv[10])	#Length of nematic pointer lines, should be about 3x the director line length
myAspect=sys.argv[11]		#'auto' - reshapes into square graph or 'equal' keeps whatever aspect ratio the true values
keepFrames=int(sys.argv[12])	#0=don't keep (delete) frames; 1=keep frames

###########################################################
### Format and style
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
myMap=ed.plasma
# cool = lsc.from_list("", [capri,ruby])
# deepsea = lsc.from_list("", [purple,ceruleandarker,limegreen])
# Adjust line width
myLW=1.0
#Animation stuff
bitrate=5000
framerate=4		#Number of frames per second in the output video
codec='libx264'		#Other options include mpeg4
suffix='.mp4'
FS=25

###########################################################
### Initialize
###########################################################
if avdim=='x':
	dim=0
	d1=1
	d2=2
elif avdim=='y':
	dim=1
	d1=0
	d2=2
elif avdim=='z':
	dim=2
	d1=0
	d2=1
else:
	print( "avdim must be 'x', 'y' or 'z' - not %s"%avdim )
	exit()

###########################################################
### Read json
###########################################################
if not os.path.isfile(inputName):
	print("%s not found."%inputName)
	exit()
with open(inputName, 'r') as f:
  input = json.load(f)
xyzSize=array([30,30,1])
dim=2
if "domain" in input:
	xyzSize[0]=input['domain'][0]
	dim=1
	if(len(input['domain'])>1):
		xyzSize[1]=input['domain'][1]
		dim=2
	else:
		xyzSize[1]=1
	if(len(input['domain'])>2):
		xyzSize[2]=input['domain'][2]
		dim=3
	else:
		xyzSize[2]=1
# Data
XYZ = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
DIR = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
XY = zeros(shape=(2,xyzSize[d1],xyzSize[d2]),dtype=float)
S = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
CHARGE = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
ANGLE = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)

### Setup the animation
# Figure
fig1 = plt.figure(1)
#Create the colorbar
CS3 = imshow(AVS.T,cmap=myMap,vmin=-0.5, vmax=0.5,aspect=myAspect)			#pcolor() sucks this is way better
cb=colorbar(CS3,shrink=float(xyzSize[d2])/float(xyzSize[d1]))
cb.ax.set_ylabel(r'Topological charge, $q$', fontsize = FS)
# Make labels
if avdim=='x':
	labX='y'
	labY='z'
elif avdim=='y':
	labX='x'
	labY='z'
elif avdim=='z':
	labX='x'
	labY='y'
#TODO: make this work for desynced avs.dat and topocharge.dat files

###########################################################
### Read the data for animation
###########################################################
print( 'Read file for figures ...' )
if not os.path.isfile(directorData):
	print("%s not found."%directorData)
	exit()
datainfile = open(directorData,"r")
if not os.path.isfile(topoData):
	print("%s not found."%topoData)
	exit()
topoinfile = open(topoData, "r")
print( '\tToss headers ...' )
for i in range(13):
	#toss header
	line = datainfile.readline()
	line2 = topoinfile.readline()
	#print line

print( '\tRead data ...' )
i=0
j=0
n=-1
while datainfile:
	i=i+1
	line = datainfile.readline()
	topoLine = topoinfile.readline()
	if( not line or not topoLine ):
		break
	else:
		if j>=start:
			#avS field
			t,Qx,Qy,Qz,Vx,Vy,Vz,s = line.split("\t",8)
			XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx) + 0.5
			XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy) + 0.5
			XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz) + 0.5
			DIR[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
			DIR[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
			DIR[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)
			S[int(Qx)][int(Qy)][int(Qz)] = float(s)

			#topo field
			t,Qx,Qy,Qz,C,angle = topoLine.split("\t", 6)
			CHARGE[int(Qx)][int(Qy)] = float(C)
			ANGLE[int(Qx)][int(Qy)] = float(angle)

	if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
		j=j+1
		if j>finish:
			break
		if j<start or j>finish:
			print( 'Toss %d'%j )
			aaa=0
		else:
			print( 'Work %d'%j )
			#Sum
			if avdim=='x':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][y][z]=MEAN[0][y][z]+DIR[0][x][y][z]
							MEAN[1][y][z]=MEAN[1][y][z]+DIR[1][x][y][z]
							MEAN[2][y][z]=MEAN[2][y][z]+DIR[2][x][y][z]
							AVS[y][z]=AVS[y][z]+S[x][y][z]
				for y in range(xyzSize[1]):
					for z in range(xyzSize[2]):
						for i in range(3):
							MEAN[i][y][z]/=xyzSize[0]
						AVS[y][z]/=xyzSize[0]
			elif avdim=='y':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][x][z]=MEAN[0][x][z]+DIR[0][x][y][z]
							MEAN[1][x][z]=MEAN[1][x][z]+DIR[1][x][y][z]
							MEAN[2][x][z]=MEAN[2][x][z]+DIR[2][x][y][z]
						AVS[x][z]=AVS[x][z]+S[x][y][z]
				for x in range(xyzSize[0]):
					for z in range(xyzSize[2]):
						for i in range(3):
							MEAN[i][x][z]/=xyzSize[1]
						AVS[x][z]/=xyzSize[1]
			elif avdim=='z':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for z in range(xyzSize[2]):
							MEAN[0][x][y]=MEAN[0][x][y]+DIR[0][x][y][z]
							MEAN[1][x][y]=MEAN[1][x][y]+DIR[1][x][y][z]
							MEAN[2][x][y]=MEAN[2][x][y]+DIR[2][x][y][z]
							AVS[x][y]=AVS[x][y]+S[x][y][z]
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						for i in range(3):
							MEAN[i][x][y]/=xyzSize[2]
						AVS[x][y]/=xyzSize[2]
			for x in range(xyzSize[d1]):
				for y in range(xyzSize[d2]):
					MAG[x][y]=sqrt( MEAN[0][x][y]**2+MEAN[1][x][y]**2+MEAN[2][x][y]**2 )

			#Save the instantaneous or current velocity field frame
			# Make Mesh
			if avdim=='x':
				for y in range(xyzSize[1]):
					for z in range(xyzSize[2]):
						XY[0][y][z]=XYZ[d1][0][y][z]
						XY[1][y][z]=XYZ[d2][0][y][z]
			elif avdim=='y':
				for x in range(xyzSize[0]):
					for z in range(xyzSize[2]):
						XY[0][x][z]=XYZ[d1][x][0][z]
						XY[1][x][z]=XYZ[d2][x][0][z]
			elif avdim=='z':
				for x in range(xyzSize[0]):
					for y in range(xyzSize[1]):
						XY[0][x][y]=XYZ[d1][x][y][0]
						XY[1][x][y]=XYZ[d2][x][y][0]

			# Save frame
			n=n+1
			fig1 = plt.figure(1)
			plt.cla()
			for x in range(xyzSize[d1]):
				for y in range(xyzSize[d2]):
					if( x%qx==0 and y%qy==0 ):
						#plot the field data
						plot( [ XY[0][x][y]-c*MEAN[d1][x][y],XY[0][x][y]+c*MEAN[d1][x][y] ],[ XY[1][x][y]-c*MEAN[d2][x][y],XY[1][x][y]+c*MEAN[d2][x][y] ],color=myMap(CHARGE[x][y]+0.5,1),linewidth=myLW )

						#plot +1/2 defects as lines 
						if (CHARGE[x][y] > 1.0e-6):
							plt.arrow(XY[0][x][y], XY[1][x][y], 0.5*c1*math.cos(ANGLE[x][y]), 0.5*c1*math.sin(ANGLE[x][y]), color='k')
							plt.arrow(XY[0][x][y], XY[1][x][y], -0.5*c1*math.cos(ANGLE[x][y]), -0.5*c1*math.sin(ANGLE[x][y]), color='k')
						#plot -1/2 defects with three lines (per symmetry)
						if (CHARGE[x][y] < -1.0e-6):
							plt.arrow(XY[0][x][y], XY[1][x][y], 0.5*c1*math.cos(ANGLE[x][y]), 0.5*c1*math.sin(ANGLE[x][y]), color='k')
							plt.arrow(XY[0][x][y], XY[1][x][y], 0.5*c1*math.cos(ANGLE[x][y] + math.pi*2/3), 0.5*c1*math.sin(ANGLE[x][y] + math.pi*2/3), color='k')
							plt.arrow(XY[0][x][y], XY[1][x][y], 0.5*c1*math.cos(ANGLE[x][y] - math.pi*2/3), 0.5*c1*math.sin(ANGLE[x][y] - math.pi*2/3), color='k')
			xlabel(r'$%s$'%labX, fontsize = FS)
			ylabel(r'$%s$'%labY, fontsize = FS)
			plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
			name='frame%04d.png'%(n)
			# savefig( name )
			savefig( name,bbox_inches='tight',pad_inches=0 )

			#Zero matrix
		DIR= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
		MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
		AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
		i=0
datainfile.close()

###########################################################
### Plot data
###########################################################
print( "Plot data ..." )

#Animate
print( "\tAnimating ..." )
name='2Dcharge_animation%s'%suffix
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
	myCommand="rm frame*.png"
	call(myCommand,shell=True)
