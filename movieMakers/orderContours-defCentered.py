from pylab import *
from numpy import ma
from subprocess import call
import sys
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.axes_grid1 import make_axes_locatable

from defectHandler import extractData

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
avdim = sys.argv[7]			# Dimension to average over
myAspect=sys.argv[8]		#'auto' - reshapes into square graph or 'equal' keeps whatever aspect ratio the true values
keepFrames=int(sys.argv[9])	#0=don't keep (delete) frames; 1=keep frames
defectData = sys.argv[10]		# Name of the defect data ("" if no defect data)
radius = int(sys.argv[11])		# Radius around defect to show snapshots for

# defect handling if needed
defects = None
print("Loading defects for rendering")
LOADDEFECTS = True

_, defects = extractData(defectData, np.array([xyzSize[0], xyzSize[1], xyzSize[2]]))
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
myLW=5.0
# adjust line length

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
			S[int(Qy)][int(Qx)][int(Qz)] = float(s) # fucking stupid hack

	if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
		j=j+1
		if j>finish:
			break
		if j<start or j>finish:
			print( 'Toss %d'%j )
			aaa=0
		else:
			print( 'Work %d'%j )

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

			# loop through defects
			for d, defect in enumerate(defects[j-1]):
				print(f"\tDefect {d+1}/{len(defects[j-1])}")
				fig1 = plt.figure(1)
				plt.cla()

				localS = np.zeros(shape=(2*radius, 2*radius), dtype=float) # S locally centered about this defect
				#assign to local region
				xDef = int(defect.pos[0])
				yDef = int(defect.pos[1])
				for x in range(-radius, radius):
					for y in range(-radius, radius):
						# apply PBC 
						globalX = np.mod(yDef + x, xyzSize[0])
						globalY = np.mod(xDef + y, xyzSize[1])
						localS[radius+x][radius+y] = S[globalX][globalY][0]

				X, Y = np.meshgrid(np.arange(xDef-radius, xDef+radius), np.arange(yDef-radius, yDef+radius))
				
				plt.contourf(X, Y, localS[:,:], origin='lower', cmap=myMap, vmin=0, vmax=1)
				defect.drawDefect(1, myLW)

				xlabel(r'$%s$'%labX, fontsize = FS)
				ylabel(r'$%s$'%labY, fontsize = FS)
				name='frame%04d-def%04d.png'%(n, d)
				# savefig( name )

				# uncomment below for snapshots
				# plt.axis('off') 
				plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

				savefig( name,bbox_inches='tight',pad_inches=0 )

			#Zero matrix
		DIR= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
		MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
		AVS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
		i=0
infile.close()
