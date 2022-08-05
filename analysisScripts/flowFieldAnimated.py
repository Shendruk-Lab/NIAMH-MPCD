"""
	Flow field rendering script.
	Uses defect handler (from https://github.com/Shendruk-Lab/MPCDDefectLoader)

	Originally from Tyler N. Shendruk
	Modified by Timofey Kozhukhov
"""

from pylab import *
from subprocess import call
import sys

from defectHandler import getDefectData

###########################################################
### Plots 2D velocity field averaging over user defined direction
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
dataName = sys.argv[1]		# path to the data (flowfield.dat)
xyzSize[0] = int(sys.argv[2])	# System size
xyzSize[1] = int(sys.argv[3])	# System size
xyzSize[2] = int(sys.argv[4])	# System size
start = int(sys.argv[5])		# Average after this number
finish = int(sys.argv[6])		# Average before this number
qx = int(sys.argv[7])		# Only show every qx arrow in x
qy = int(sys.argv[8])		# Only show every qy arrow in y
avdim = sys.argv[9]			# Dimension to average over
myAspect=sys.argv[10]		#'auto' - reshapes into square graph or 'equal' keeps whatever aspect ratio the true values
keepFrames=int(sys.argv[11])	#0=don't keep (delete) frames; 1=keep frames
savePDF=int(sys.argv[12]) # 1 for saving transparent pdfs (for papers), 0 for none
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
myMap=ed.deepsea

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
VEL = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
currentMEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
currentMAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
XY = zeros(shape=(2,xyzSize[d1],xyzSize[d2]),dtype=float)

###########################################################
### Read the data for min/max
###########################################################
print( 'Read file for min/max ...' )
file = dataName
infile = open(file,"r")
minV=99999999999999.0
maxV=0.0

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
  if( len(line)!= 57):
    break
  else:
    Qx,Qy,Qz,Vx,Vy,Vz = line.split("\t",6)
    XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx) + 0.5
    XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy) + 0.5
    XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz) + 0.5
    VEL[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
    VEL[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
    VEL[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)

  if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
    j=j+1
    print( 'Work %d'%j )

    if j>finish:
      break
    if j<start or j>finish:
      print( 'Toss %d'%j )
    else:
      #print( 'Work %d'%j )
      #Sum
      if avdim=='x':
        for x in range(xyzSize[0]):
          for y in range(xyzSize[1]):
            for z in range(xyzSize[2]):
              MEAN[0][y][z]=MEAN[0][y][z]+VEL[0][x][y][z]
              MEAN[1][y][z]=MEAN[1][y][z]+VEL[1][x][y][z]
              MEAN[2][y][z]=MEAN[2][y][z]+VEL[2][x][y][z]
              currentMEAN[0][y][z]=currentMEAN[0][y][z]+VEL[0][x][y][z]
              currentMEAN[1][y][z]=currentMEAN[1][y][z]+VEL[1][x][y][z]
              currentMEAN[2][y][z]=currentMEAN[2][y][z]+VEL[2][x][y][z]
        for y in range(xyzSize[1]):
          for z in range(xyzSize[2]):
            for i in range(3):
              currentMEAN[i][y][z]/=xyzSize[0]
      elif avdim=='y':
        for x in range(xyzSize[0]):
          for y in range(xyzSize[1]):
            for z in range(xyzSize[2]):
              MEAN[0][x][z]=MEAN[0][x][z]+VEL[0][x][y][z]
              MEAN[1][x][z]=MEAN[1][x][z]+VEL[1][x][y][z]
              MEAN[2][x][z]=MEAN[2][x][z]+VEL[2][x][y][z]
              currentMEAN[0][x][z]=currentMEAN[0][x][z]+VEL[0][x][y][z]
              currentMEAN[1][x][z]=currentMEAN[1][x][z]+VEL[1][x][y][z]
              currentMEAN[2][x][z]=currentMEAN[2][x][z]+VEL[2][x][y][z]
        for x in range(xyzSize[0]):
          for z in range(xyzSize[2]):
            for i in range(3):
              currentMEAN[i][x][z]/=xyzSize[1]
      elif avdim=='z':
        for x in range(xyzSize[0]):
          for y in range(xyzSize[1]):
            for z in range(xyzSize[2]):
              MEAN[0][x][y]=MEAN[0][x][y]+VEL[0][x][y][z]
              MEAN[1][x][y]=MEAN[1][x][y]+VEL[1][x][y][z]
              MEAN[2][x][y]=MEAN[2][x][y]+VEL[2][x][y][z]
              currentMEAN[0][x][y]=currentMEAN[0][x][y]+VEL[0][x][y][z]
              currentMEAN[1][x][y]=currentMEAN[1][x][y]+VEL[1][x][y][z]
              currentMEAN[2][x][y]=currentMEAN[2][x][y]+VEL[2][x][y][z]
        for x in range(xyzSize[0]):
          for y in range(xyzSize[1]):
            for i in range(3):
              currentMEAN[i][x][y]/=xyzSize[2]
        for x in range(xyzSize[d1]):
          for y in range(xyzSize[d2]):
            currentMAG[x][y]=sqrt( currentMEAN[0][x][y]**2+currentMEAN[1][x][y]**2+currentMEAN[2][x][y]**2 )

        for x in range(xyzSize[d1]):
          for y in range(xyzSize[d2]):
            if currentMAG[x][y]>maxV:
              maxV=currentMAG[x][y]
            elif currentMAG[x][y]<minV:
              minV=currentMAG[x][y]

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

    ###########################################################
    ### Plot the frame
    ###########################################################
    # Save frame
    n=n+1
    fig1, axes = plt.subplots(nrows=1, ncols=1)
    #plt.cla()

    if j==1:
      #Setup the density image
      plt.subplot(1,1,1)
      plt.cla()
      #Setup the velocity image
      quiv = quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], currentMEAN[d1][::qx, ::qy], currentMEAN[d2][::qx, ::qy] )
      velImage = imshow(currentMAG.T,cmap=myMap,origin='lower',aspect=myAspect,vmin=minV,vmax=maxV)
      velCB = colorbar()
      velCB.ax.set_ylabel(r'Velocity, $\left|\vec{v}\right|$', fontsize = FS)
      xlabel(r'$%s$'%labX, fontsize = FS)
      ylabel(r'$%s$'%labY, fontsize = FS)
    else:
      #Velocity image
      plt.subplot(1,1,1)
      plt.cla()
      quiv = quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], currentMEAN[d1][::qx, ::qy], currentMEAN[d2][::qx, ::qy] )
      velImage = imshow(currentMAG.T,cmap=myMap,origin='lower',aspect=myAspect,vmin=minV,vmax=maxV)
      velCB = colorbar()
      velCB.ax.set_ylabel(r'Velocity, $\left|\vec{v}\right|$', fontsize = FS)
      fig1.canvas.draw()
      xlabel(r'$%s$'%labX, fontsize = FS)
      ylabel(r'$%s$'%labY, fontsize = FS)
    
    plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)

    # load defects and draw them as necessary
    # FIXME: only works for 2d for now, doesnt take into account d1 or d2
    if LOADDEFECTS and (j < len(defects)):
      print(f"\tDrawing defects {j}/{len(defects)-1}")
      for defect in defects[j-1]: # j is not 0 indexed reeeeee
        defect.drawDefect(0.5*(qx+qy), 2)

    name='frame%04d.png'%(n)
    namepdf='frame%04d.pdf'%(n)

    ## uncomment below for snapshots!
    plt.axis('off') 
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    plt.savefig(name, bbox_inches='tight')
    # save trans pdf
    if savePDF: plt.savefig(namepdf, transparent=True, bbox_inches='tight')

    #Zero matrix
    VEL= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
    currentMEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
    i=0
infile.close()

###########################################################
### Average the data
###########################################################
#This x and y aren't necessarily the actual x and y
print( "Average data over %d instances ..."%(j-start) )
norm=float(xyzSize[dim])
for x in range(xyzSize[d1]):
  for y in range(xyzSize[d2]):
    MEAN[0][x][y]=MEAN[0][x][y]/norm/float(j-start)
    MEAN[1][x][y]=MEAN[1][x][y]/norm/float(j-start)
    MEAN[2][x][y]=MEAN[2][x][y]/norm/float(j-start)
    MAG[x][y]=sqrt( MEAN[0][x][y]**2+MEAN[1][x][y]**2+MEAN[2][x][y]**2 )

###########################################################
### Plot data
###########################################################
print( "Plot data ..." )

# Animate
print( "\tAnimating ..." )
name='2Dvelocity_animation%s'%suffix
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
    myCommand="rm frame*.png"
    call(myCommand,shell=True)

print( "\tPlotting ..." )
fig, ax = plt.subplots()
imshow(MAG.T,cmap=myMap,origin='lower',aspect='auto')
cb=colorbar()
cb.ax.set_ylabel(r'$\left|\vec{u}\right|$', fontsize = FS)
quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], MEAN[d1][::qx, ::qy], MEAN[d2][::qx, ::qy] )
xlabel(r'$%s$'%labX, fontsize = FS)
ylabel(r'$%s$'%labY, fontsize = FS)
ax.tick_params(axis='both', which='major', labelsize=TLS)
plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
name='2D_av_%s.pdf'%avdim
print( "\t%s"%name )
savefig( name )

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(XY[0],XY[1],MAG, rstride=1,cstride=1,cmap=myMap,linewidth=0,antialiased=False,alpha=0.6)
ax.plot_wireframe(XY[0],XY[1],MAG, rstride=1,cstride=1,linewidth=0.25,color='k',alpha=1.0)
ax.set_xlabel(r'$%s$'%labX, fontsize = FS)
ax.set_ylabel(r'$%s$'%labY, fontsize = FS)
ax.set_zlabel(r'$\left|\vec{u}\right|$', fontsize = FS)
ax.tick_params(axis='both', which='major', labelsize=TLS)
name='2Dcontour_av_%s.pdf'%avdim
print( "\t%s"%name )
savefig( name )

#plt.show()
