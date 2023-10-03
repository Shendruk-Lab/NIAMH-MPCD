"""
	Orientation field (director) rendering script.
	Uses defect handler (from https://github.com/Shendruk-Lab/MPCDDefectLoader)

	Originally from Tyler N. Shendruk
	Modified by Timofey Kozhukhov
"""

from pylab import *
from subprocess import call
import sys

from defectHandler import getDefectData

###########################################################
### Plots 2D averaging over user defined direction
### Input files: directorfield.dat (OrderParaAndDirField)
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

dataName = sys.argv[1]		# Name of the data
xyzSize[0] = int(sys.argv[2])	# System size
xyzSize[1] = int(sys.argv[3])	# System size
xyzSize[2] = int(sys.argv[4])	# System size
start = int(sys.argv[5])		# Average after this number
finish = int(sys.argv[6])		# Average before this number
qx = int(sys.argv[7])		# Only show every qx arrow in x
qy = int(sys.argv[8])		# Only show every qy arrow in y
avdim = sys.argv[9]			# Dimension to average over
c = float(sys.argv[10])		#Length of director lines approx 0.5
myAspect=sys.argv[11]		#'auto' - reshapes into square graph or 'equal' keeps whatever aspect ratio the true values
keepFrames=int(sys.argv[12])	#0=don't keep (delete) frames; 1=keep frames
defectData = ""
try:
	defectData = sys.argv[13]		# Name of the defect data ("" if no defect data)
except:
	print("No defect data found")

# defect handling if needed
LOADDEFECTS = False
defects = []
if defectData != "":
	print("Loading defects for rendering")
	LOADDEFECTS = True

	defContainer = getDefectData(defectData, np.array([xyzSize[0], xyzSize[1], xyzSize[2]]))
	for defList in defContainer:
		defects.append(defList.defectList)
	print("Finished loading defects")

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
# adjust line length
c *= (qx+qy)/2

#Animation stuff
bitrate=5000
framerate=12		#Number of frames per second in the output video
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

# Figure
fig1 = plt.figure(1)
#Create the colorbar
CS3 = imshow(AVS.T,cmap=myMap,vmin=0, vmax=1,aspect=myAspect)			#pcolor() sucks this is way better
cb=colorbar(CS3,shrink=float(xyzSize[d2])/float(xyzSize[d1]))
cb.ax.set_ylabel(r'$S$', fontsize = FS)

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
file = dataName
infile = open(file,"r")
print( '\tToss header ...' )
for i in range(13):
	#toss header
	line = infile.readline()
	#print line

print( '\tRead data ...' )
i=0
j=0
n=-1
while infile:
	i=i+1
	line = infile.readline()
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

			quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], 
					c*MEAN[0][::qx, ::qy], c*MEAN[1][::qx, ::qy], 
					AVS[::qx, ::qy], cmap=myMap, clim=(0, 1), scale=50/c,
					width=0.005*myLW, headlength=0, headwidth=0, 
					headaxislength=0, pivot='middle')
					# linewidth=myLW, headaxislength=1)

			# old method - draw a bunch of individual lines. slow and painful
			# for x in range(xyzSize[d1]):
			# 	for y in range(xyzSize[d2]):
			# 		if( x%qx==0 and y%qy==0 ):
			# 			plot( [ XY[0][x][y]-c*MEAN[d1][x][y],XY[0][x][y]+c*MEAN[d1][x][y] ],[ XY[1][x][y]-c*MEAN[d2][x][y],XY[1][x][y]+c*MEAN[d2][x][y] ],color=myMap(AVS[x][y],1),linewidth=myLW )

			# load defects and draw them as necessary
			# FIXME: only works for 2d for now, doesnt take into account d1 or d2
			if LOADDEFECTS and (j < len(defects)):
				print(f"\tDrawing defects {j}/{len(defects)-1}")
				for defect in defects[j-1]: # j is not 0 indexed reeeeee
					defect.drawDefect(c*2, myLW*2)

			xlabel(r'$%s$'%labX, fontsize = FS)
			ylabel(r'$%s$'%labY, fontsize = FS)
			plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
			name='frame%04d.png'%(n)
			# savefig( name )

			# uncomment below for snapshots
			plt.axis('off') 
			plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

			savefig( name,bbox_inches='tight',pad_inches=0 )

			#Zero matrix
		DIR= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
		MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
		AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
		i=0
infile.close()

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
