"""
	Flow field rendering script.
	Uses defect handler (from https://github.com/Shendruk-Lab/MPCDDefectLoader)

	Originally from Tyler N. Shendruk
	Modified by Timofey Kozhukhov
"""

from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from subprocess import call
import sys
import os
import json
import argparse

from defectHandler import getDefectData

###########################################################
### Set up argsparse
###########################################################
parser = argparse.ArgumentParser(description='Flow field rendering script.')
parser.add_argument("dataname", type=str, help="Path to the data (should be "
                                               "flowfield.dat)")
parser.add_argument("inputname", type=str, help="Path to input .json file")
parser.add_argument("start", type=int, help="Starting timestep for averaging")
parser.add_argument("finish", type=int, help="Finishing timestep for averaging")
parser.add_argument("--qx", type=int, help="Only show every qx arrow in x",
                    default=1)
parser.add_argument("--qy", type=int, help="Only show every qy arrow in y",
                    default=1)
parser.add_argument("avdim", type=str, help="Dimension to average over")
parser.add_argument("-a", "--myAspect", type=str, help="'auto' or 'equal'",
                    default="auto")
parser.add_argument("-k", "--keepFrames", type=int,
                    help="0=don't keep (delete) frames; 1=keep frames",
                    default=0)
parser.add_argument("-p", "--savePDF", type=int,
                    help="1 saves transparent pdfs for papers, 0 for none",
                    default=0)
parser.add_argument("-d", "--defectData", type=str,
                    help="Path to defect data (if any)", default="")
args = parser.parse_args()

###########################################################
### Read arguments
###########################################################
print("Arguments:")
for arg, value in vars(args).items():
	print(f"\t{arg}: {value}")
dataName = args.dataname
inputName = args.inputname
start = args.start
finish = args.finish
qx = args.qx
qy = args.qy
avdim = args.avdim
myAspect = args.myAspect
keepFrames = args.keepFrames
savePDF = args.savePDF
defectData = args.defectData

###########################################################
### Format and style
###########################################################
# Use our custom style and colours
plt.style.use('shendrukGroupStyle')
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
VEL = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
currentMEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
currentMAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
XY = zeros(shape=(2,xyzSize[d1],xyzSize[d2]),dtype=float)

###########################################################
### defect handling if needed
###########################################################
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
### Read the data for min/max
###########################################################
print( 'Reading file for min/max ...' )
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
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
infile = open(file,"r")
# Toss header
for i in range(13):
  line = infile.readline()

i=0
j=0
n=-1
while infile:
  i=i+1
  line = infile.readline()
  if not line:
    break
  else:
    split = line.split("\t")
    if len(split) == 6:
      Qx,Qy,Qz,Vx,Vy,Vz = split
    elif len(split) == 7:
      _,Qx,Qy,Qz,Vx,Vy,Vz = split # to take into account old vs new format
    else:
      break # error handling
    XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx) + 0.5
    XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy) + 0.5
    XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz) + 0.5
    VEL[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
    VEL[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
    VEL[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)

  if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
    j=j+1
    if j>finish:
      break
    if j<start or j>finish:
      pass
    else:
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
    if j==1:
      #Setup the density image
      plt.subplot(1,1,1)
      plt.cla()
      #Setup the velocity image
      quiv = quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], currentMEAN[d1][::qx, ::qy], currentMEAN[d2][::qx, ::qy] )
      velImage = imshow(currentMAG.T,cmap=myMap,origin='lower',aspect=myAspect,vmin=minV,vmax=maxV)
      velCB = colorbar()
      velCB.ax.set_ylabel(r'Velocity, $\left|\vec{v}\right|$')
      xlabel(r'$%s$'%labX)
      ylabel(r'$%s$'%labY)
    else:
      #Velocity image
      plt.subplot(1,1,1)
      plt.cla()
      quiv = quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], currentMEAN[d1][::qx, ::qy], currentMEAN[d2][::qx, ::qy] )
      velImage = imshow(currentMAG.T,cmap=myMap,origin='lower',aspect=myAspect,vmin=minV,vmax=maxV)
      velCB = colorbar()
      velCB.ax.set_ylabel(r'Velocity, $\left|\vec{v}\right|$')
      fig1.canvas.draw()
      xlabel(r'$%s$'%labX)
      ylabel(r'$%s$'%labY)    
    plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)

    # load defects and draw them as necessary
    # FIXME: only works for 2d for now, doesnt take into account d1 or d2
    if LOADDEFECTS and (j < len(defects)):
      print(f"Drawing defects {j}/{len(defects)-1}")
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
print( "Averaging data over %d instances ..."%(j-start) )
norm=float(xyzSize[dim])
for x in range(xyzSize[d1]):
  for y in range(xyzSize[d2]):
    MEAN[0][x][y]=MEAN[0][x][y]/norm/float(j-start)
    MEAN[1][x][y]=MEAN[1][x][y]/norm/float(j-start)
    MEAN[2][x][y]=MEAN[2][x][y]/norm/float(j-start)
    MAG[x][y]=sqrt( MEAN[0][x][y]**2+MEAN[1][x][y]**2+MEAN[2][x][y]**2 )

###########################################################
### Visualize data
###########################################################
# Animate
print( "Animating ..." )
name='2Dvelocity_animation%s'%suffix
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
    myCommand="rm frame*.png"
    call(myCommand,shell=True)

print( "Plotting ..." )
fig, ax = plt.subplots()
imshow(MAG.T,cmap=myMap,origin='lower',aspect='auto')
cb=colorbar()
cb.ax.set_ylabel(r'$\left|\vec{u}\right|$')
quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], MEAN[d1][::qx, ::qy], MEAN[d2][::qx, ::qy] )
xlabel(r'$%s$'%labX)
ylabel(r'$%s$'%labY)
ax.tick_params(axis='both', which='major')
plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
name='2D_av_%s.pdf'%avdim
savefig( name )

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(XY[0],XY[1],MAG, rstride=1,cstride=1,cmap=myMap,linewidth=0,antialiased=False,alpha=0.6)
ax.plot_wireframe(XY[0],XY[1],MAG, rstride=1,cstride=1,linewidth=0.25,color='k',alpha=1.0)
ax.set_xlabel(r'$%s$'%labX)
ax.set_ylabel(r'$%s$'%labY)
ax.set_zlabel(r'$\left|\vec{u}\right|$')
ax.tick_params(axis='both', which='major')
name='2Dcontour_av_%s.pdf'%avdim
savefig( name )

#plt.show()
