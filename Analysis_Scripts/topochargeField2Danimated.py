from pylab import *
from numpy import ma
from subprocess import call
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.axes_grid1 import make_axes_locatable

###########################################################
### Plots 2D averaging over user defined direction
###########################################################

###########################################################
### Read arguments
###########################################################
FS =25
TLS = 20		# Tick label size
xyzSize=zeros( 3,dtype=int )
print( "Arguments:" )
for arg in sys.argv:
	print( "\t" + arg )
directorData = sys.argv[1]		# Name of the field data
topoData = sys.argv[2]			# Name of the topo data
xyzSize[0] = int(sys.argv[3])	# System size
xyzSize[1] = int(sys.argv[4])	# System size
xyzSize[2] = int(sys.argv[5])	# System size
start = int(sys.argv[6])		# Average after this number
finish = int(sys.argv[7])		# Average before this number
qx = int(sys.argv[8])		# Only show every qx arrow in x
qy = int(sys.argv[9])		# Only show every qy arrow in y
avdim = sys.argv[10]			# Dimension to average over
c = float(sys.argv[11])		#Length of director lines approx 0.5
myAspect=sys.argv[12]		#'auto' - reshapes into square graph or 'equal' keeps whatever aspect ratio the true values
keepFrames=int(sys.argv[13])	#0=don't keep (delete) frames; 1=keep frames

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

# Data
XYZ = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
DIR = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
XY = zeros(shape=(2,xyzSize[d1],xyzSize[d2]),dtype=float)
S = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
CHARGE = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)

# Figure
fig1 = plt.figure(1)
#Create the colorbar
CS3 = imshow(AVS.T,cmap=myMap,vmin=-0.5, vmax=0.5,aspect=myAspect)			#pcolor() sucks this is way better
cb=colorbar(CS3,shrink=float(xyzSize[d2])/float(xyzSize[d1]))
cb.ax.set_ylabel(r'$S$', fontsize = FS)

#TODO: make this work for desynced avs.dat and topocharge.dat files

###########################################################
### Read the data for animation
###########################################################

### Setup the animation
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

print( 'Read file for figures ...' )
file = directorData
datainfile = open(file,"r")
file = topoData
topoinfile = open(file, "r")
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
			t,Qx,Qy,Qz,C = topoLine.split("\t", 5)
			CHARGE[int(Qx)][int(Qy)] = float(C)

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
						plot( [ XY[0][x][y]-c*MEAN[d1][x][y],XY[0][x][y]+c*MEAN[d1][x][y] ],[ XY[1][x][y]-c*MEAN[d2][x][y],XY[1][x][y]+c*MEAN[d2][x][y] ],color=myMap(CHARGE[x][y]+0.5,1),linewidth=myLW )
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
name='2Dorientation_animation%s'%suffix
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
	myCommand="rm frame*.png"
	call(myCommand,shell=True)
