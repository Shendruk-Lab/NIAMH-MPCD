###########################################################
### Animates swimmers and 2D fields, averaging over user defined direction
###########################################################

###########################################################
### Imports
###########################################################
from pylab import *
from subprocess import call
import os
import json

###########################################################
### Read arguments
###########################################################
print( "Arguments:" )
for arg in sys.argv:
  print( "\t" + arg )
dataPath = sys.argv[1]		  # Name of the data
inputName = sys.argv[2]			# Input json file to read inputs
start = int(sys.argv[3])		# Average after this number
finish = int(sys.argv[4])		# Average before this number
qx = int(sys.argv[5])   		# Only show every qx arrow in x
qy = int(sys.argv[6])		    # Only show every qy arrow in y
avdim = sys.argv[7]			    # Dimension to average over
fieldType = sys.argv[8]		  # fieldType can be "vel", "vort" or "nem" for velocity, vorticity or director field, respectively
myAspect=sys.argv[9]    		#'auto' - reshapes into square graph or 'equal' keeps whatever aspect ratio the true values
keepFrames=int(sys.argv[10])	#0=don't keep (delete) frames; 1=keep frames

###########################################################
### Style/formating stuff
###########################################################
# Use our custom style and colours
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed
# Colours
swimmerMap=ed.viridis
if(fieldType=="vel"):
  fieldMap = ed.truncate_colormap(ed.viridis, minval=0.1, maxval=0.7)
elif(fieldType=="vor"):
  fieldMap = ed.bombpops
elif(fieldType=="nem"):
  fieldMap = ed.plasma
else:
  print("Field type %s not known."%(fieldType))
  exit()
#Animation
bitrate=5000
framerate=20		#Number of frames per second in the output video
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

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
numSw=0
if "nSwim" in input:
    numSw=input['nSwim']
dipole=1
if "dsSwim" in input:
    dipole=input['dsSwim']

###########################################################
### Set arguments
###########################################################
tailFreq=1.0/3.5
tailWaveLength=(4.0/10.0)*2.0*(1.0+dipole)
hidgeonLength=1.0
tailRad=1.0
fieldType=fieldType.lower()
if(fieldType=="v" or fieldType=="vel" or fieldType=="velocity"):
  fieldType="vel"
elif(fieldType=="w" or fieldType=="vor" or fieldType=="vort" or fieldType=="vorticity"):
  fieldType="vor"
elif(fieldType=="n" or fieldType=="nem" or fieldType=="nematic" or fieldType=="dir" or fieldType=="director"):
  fieldType="nem"

###########################################################
### Initialize
###########################################################
velFieldName="%s/flowfield.dat"%dataPath
dirFieldName="%s/directorfield.dat"%dataPath
swimmerName="%s/swimmers.dat"%dataPath
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
VORTZ = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
XY = zeros(shape=(2,xyzSize[d1],xyzSize[d2]),dtype=float)
if(fieldType=="nem"):
  XYZorder = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
  XYorder = zeros(shape=(2,xyzSize[d1],xyzSize[d2]),dtype=float)
  DIR = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
  currentDIR = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
  S = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
  currentS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
  dirL=0.25*sqrt(qx*qx+qy*qy)

# Setup figure
fig1, axes = plt.subplots(nrows=1, ncols=1)

###########################################################
### Read the data for min/max
###########################################################
print( 'Reading data ...' )
if not os.path.isfile(swimmerName):
	print("%s not found."%swimmerName)
	exit()
swimFile = open(swimmerName,"r")
if not os.path.isfile(velFieldName):
	print("%s not found."%velFieldName)
	exit()
velFile = open(velFieldName,"r")
minV=99999999999999.0
maxV=0.0
minW=0.0
maxW=0.0
if(fieldType=="nem"):
  if not os.path.isfile(dirFieldName):
    print("%s not found."%dirFieldName)
    exit()
  nemFile = open(dirFieldName,"r")

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
xyzF=zeros( 3,dtype=float )
for i in range(3):
  xyzF[i]=float(xyzSize[i])

# Toss headers
for i in range(13):
  line = velFile.readline()
for i in range(13):
  line = swimFile.readline()
if(fieldType=="nem"):
  for i in range(13):
    line = nemFile.readline()
# Read data
i=0
j=0
n=-1
while( velFile ):
  i=i+1
  ###########################################################
  ### Read the order
  ###########################################################
  if(fieldType=="nem"):
    line = nemFile.readline()
    if( not line ):
      break
    else:
      if j>=start:
        t,Qx,Qy,Qz,Vx,Vy,Vz,s = line.split("\t",8)
        XYZorder[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx) + 0.5
        XYZorder[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy) + 0.5
        XYZorder[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz) + 0.5
        DIR[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
        DIR[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
        DIR[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)
        S[int(Qx)][int(Qy)][int(Qz)] = float(s)

  ###########################################################
  ### Read the velocity
  ###########################################################
  line = velFile.readline()
  if( len(line)!= 70):
    break
  else:
    time,Qx,Qy,Qz,Vx,Vy,Vz = line.split("\t",7)
    XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx) + 0.5
    XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy) + 0.5
    XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz) + 0.5
    VEL[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
    VEL[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
    VEL[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)

  if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
    ###########################################################
    ### Average the fields
    ###########################################################
    j=j+1
    if j>finish:
      break
    if j<start or j>finish:
      pass
    else:
      ###########################################################
      ### Order
      ###########################################################
      if(fieldType=="nem"):
        if avdim=='x':
          for x in range(xyzSize[0]):
            for y in range(xyzSize[1]):
              for z in range(xyzSize[2]):
                currentDIR[0][y][z]=currentDIR[0][y][z]+DIR[0][x][y][z]
                currentDIR[1][y][z]=currentDIR[1][y][z]+DIR[1][x][y][z]
                currentDIR[2][y][z]=currentDIR[2][y][z]+DIR[2][x][y][z]
                currentS[y][z]=currentS[y][z]+S[x][y][z]
          for y in range(xyzSize[1]):
            for z in range(xyzSize[2]):
              for d in range(3):
                currentDIR[d][y][z]/=xyzSize[0]
                currentS[y][z]/=xyzSize[0]
        elif avdim=='y':
          for x in range(xyzSize[0]):
            for y in range(xyzSize[1]):
              for z in range(xyzSize[2]):
                currentDIR[0][x][z]=currentDIR[0][x][z]+DIR[0][x][y][z]
                currentDIR[1][x][z]=currentDIR[1][x][z]+DIR[1][x][y][z]
                currentDIR[2][x][z]=currentDIR[2][x][z]+DIR[2][x][y][z]
                currentS[x][z]=currentS[x][z]+S[x][y][z]
          for x in range(xyzSize[0]):
            for z in range(xyzSize[2]):
              for d in range(3):
                currentDIR[d][x][z]/=xyzSize[1]
                currentS[x][z]/=xyzSize[1]
        elif avdim=='z':
          for x in range(xyzSize[0]):
            for y in range(xyzSize[1]):
              for z in range(xyzSize[2]):
                currentDIR[0][x][y]=currentDIR[0][x][y]+DIR[0][x][y][z]
                currentDIR[1][x][y]=currentDIR[1][x][y]+DIR[1][x][y][z]
                currentDIR[2][x][y]=currentDIR[2][x][y]+DIR[2][x][y][z]
                currentS[x][y]=currentS[x][y]+S[x][y][z]
          for x in range(xyzSize[0]):
            for y in range(xyzSize[1]):
              for d in range(3):
                currentDIR[d][x][y]/=xyzSize[2]
                currentS[x][y]/=xyzSize[2]
        # Make Mesh
        if avdim=='x':
          for y in range(xyzSize[1]):
            for z in range(xyzSize[2]):
              XYorder[0][y][z]=XYZorder[d1][0][y][z]
              XYorder[1][y][z]=XYZorder[d2][0][y][z]
        elif avdim=='y':
          for x in range(xyzSize[0]):
            for z in range(xyzSize[2]):
              XYorder[0][x][z]=XYZorder[d1][x][0][z]
              XYorder[1][x][z]=XYZorder[d2][x][0][z]
        elif avdim=='z':
          for x in range(xyzSize[0]):
            for y in range(xyzSize[1]):
              XYorder[0][x][y]=XYZorder[d1][x][y][0]
              XYorder[1][x][y]=XYZorder[d2][x][y][0]

      ###########################################################
      ### Velocity
      ###########################################################
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
            for d in range(3):
              currentMEAN[d][y][z]/=xyzSize[0]
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
            for d in range(3):
              currentMEAN[d][x][z]/=xyzSize[1]
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
            for d in range(3):
              currentMEAN[d][x][y]/=xyzSize[2]
      # Calculate magnitude
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
      ### Vorticity
      ###########################################################
      if(fieldType=="vor"):
        #This x and y aren't necessarily the actual x and y
        dx=1.0
        dy=1.0
        #Corners both x and y are forward or reverse derivatives
        #0,0
        vydx=(currentMEAN[d2][1][0]-currentMEAN[d2][0][0])/dx
        vxdy=(currentMEAN[d1][0][1]-currentMEAN[d1][0][0])/dy
        VORTZ[0][0] =vydx-vxdy
        #N,0
        vydx=(currentMEAN[d2][-1][0]-currentMEAN[d2][-2][0])/dx
        vxdy=(currentMEAN[d1][-1][0]-currentMEAN[d1][-1][1])/dy
        VORTZ[-1][0] =vydx-vxdy
        #N,N
        vydx=(currentMEAN[d2][-1][-1]-currentMEAN[d2][-2][-1])/dx
        vxdy=(currentMEAN[d1][-1][-1]-currentMEAN[d1][-1][-2])/dy
        VORTZ[-1][-1] =vydx-vxdy
        #0,N
        vydx=(currentMEAN[d2][1][-1]-currentMEAN[d2][0][-1])/dx
        vxdy=(currentMEAN[d1][0][-1]-currentMEAN[d1][0][-2])/dy
        VORTZ[-1][-1] =vydx-vxdy
        #Edges are a forward or reverse derivative with centred derivatives
        for x in range(1,xyzSize[d1]-1):
          vydx=0.5*(currentMEAN[d2][x+1][0]-currentMEAN[d2][x-1][0])/dx
          vxdy=(currentMEAN[d1][x][1]-currentMEAN[d1][x][0])/dy
          VORTZ[x][0] =vydx-vxdy
          vydx=0.5*(currentMEAN[d2][x+1][-1]-currentMEAN[d2][x-1][-1])/dx
          vxdy=(currentMEAN[d1][x][-1]-currentMEAN[d1][x][-2])/dy
          VORTZ[x][-1] =vydx-vxdy
        for y in range(1,xyzSize[d2]-1):
          vydx=(currentMEAN[d2][1][y]-currentMEAN[d2][0][y])/dx
          vxdy=0.5*(currentMEAN[d1][0][y+1]-currentMEAN[d1][0][y-1])/dy
          VORTZ[0][y] =vydx-vxdy
          vydx=(currentMEAN[d2][-1][y]-currentMEAN[d2][-2][y])/dx
          vxdy=0.5*(currentMEAN[d1][-1][y+1]-currentMEAN[d1][-1][y-1])/dy
          VORTZ[-1][y] =vydx-vxdy
        #The bulk is all centred derivatives
        for x in range(1,xyzSize[d1]-1):
          for y in range(1,xyzSize[d2]-1):
            vydx=0.5*(currentMEAN[d2][x+1][y]-currentMEAN[d2][x-1][y])/dx
            vxdy=0.5*(currentMEAN[d1][x][y+1]-currentMEAN[d1][x][y-1])/dy
            VORTZ[x][y] =vydx-vxdy
        for x in range(xyzSize[d1]):
          for y in range(xyzSize[d2]):
            if VORTZ[x][y]>maxW:
              maxW=VORTZ[x][y]
            elif VORTZ[x][y]<minW:
              minW=VORTZ[x][y]

    ###########################################################
    ### Read the swimmers' positions
    ###########################################################
    H=[]
    B=[]
    T=[]
    for ns in range(numSw):
      line = swimFile.readline()
      if( not line ):
        break
      else:
        t,Hx,Hy,Hz,Hvx,Hvy,Hvz,Bx,By,Bz,Bvx,Bvy,Bvz,phase = line.split( )
        H.append( [float(Hx),float(Hy),float(Hz)] )
        B.append( [float(Bx),float(By),float(Bz)] )
        T.append( [(1.0+dipole)*B[-1][0]-dipole*H[-1][0],(1.0+dipole)*B[-1][1]-dipole*H[-1][1],(1.0+dipole)*B[-1][2]-dipole*H[-1][2]] )
    #Toss the line break
    line = swimFile.readline()
    if( not line ):
      break

    ###########################################################
    ### Plot the frame
    ###########################################################
    # Save frame
    n=n+1
    plt.subplot(1,1,1)
    #Setup the velocity image
    quiv = quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], currentMEAN[d1][::qx, ::qy], currentMEAN[d2][::qx, ::qy] )
    if(fieldType=="vel"):
      image = imshow(currentMAG.T,cmap=fieldMap,origin='lower',aspect=myAspect,vmin=minV,vmax=maxV)
      CB = colorbar()
      CB.ax.set_ylabel(r'Velocity, $\left|\vec{v}\right|$')
    elif(fieldType=="vor"):
      image = imshow(VORTZ.T,cmap=fieldMap,origin='lower',aspect=myAspect,vmin=minW,vmax=maxW)
      CB = colorbar()
      CB.ax.set_ylabel(r'Vorticity, $\omega_{%s}$'%(avdim))
    elif(fieldType=="nem"):
      for x in range(xyzSize[d1]):
        for y in range(xyzSize[d2]):
          if( x%qx==0 and y%qy==0 ):
            plot( [ XYorder[0][x][y]-dirL*currentDIR[d1][x][y],XYorder[0][x][y]+dirL*currentDIR[d1][x][y] ],[ XYorder[1][x][y]-dirL*currentDIR[d2][x][y],XYorder[1][x][y]+dirL*currentDIR[d2][x][y] ],color=fieldMap(currentS[x][y],1) )
    # Plot the swimmers
    for ns in range(numSw):
      if( fabs(H[ns][d1]-B[ns][d1])<0.5*xyzSize[d1] and fabs(H[ns][d2]-B[ns][d2])<0.5*xyzSize[d2] ):
        plot( [H[ns][d1],B[ns][d1]],[H[ns][d2],B[ns][d2]],'-',color='k',linewidth=2 )
      if( fabs(T[ns][d1]-B[ns][d1])<0.5*xyzSize[d1] and fabs(T[ns][d2]-B[ns][d2])<0.5*xyzSize[d2] ):
        plot( [T[ns][d1],B[ns][d1]],[T[ns][d2],B[ns][d2]],'-',color=swimmerMap(B[ns][dim]/xyzF[dim]),linewidth=2 )
      waveLength=(4.0/5.0)*2.0*(1.0+dipole)
      theta=arctan2(H[ns][d2]-B[ns][d2],H[ns][d1]-B[ns][d1])
      X=linspace(0,2.0*(1.0+dipole),20)
      Y=tailRad*(1.0-exp(-pow(X/hidgeonLength,2)))*sin(2.0*pi*X/tailWaveLength + float(n)*tailFreq)
      tailX=cos(theta)*X - sin(theta)*Y
      tailY=sin(theta)*X + cos(theta)*Y
      plot( B[ns][d1]-tailX,B[ns][d2]-tailY,'-',color=swimmerMap(B[ns][dim]/xyzF[dim]),linewidth=2 )
      plot( H[ns][d1],H[ns][d2],'o',color='k',fillstyle='full' )
      plot( B[ns][d1],B[ns][d2],'h',color=swimmerMap(B[ns][dim]/xyzF[dim]),fillstyle='full' )
    xlabel(r'$%s$'%labX)
    ylabel(r'$%s$'%labY)
    plt.axis(xmax=xyzSize[0], xmin=0, ymax=xyzSize[1], ymin=0)
    name='frame%04d.png'%(n)
    savefig( name )
    plt.clf()
    #Zero matrices
    VEL= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
    currentMEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
    VORTZ = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
    S = zeros( (xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
    currentS = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
    currentDIR= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
    i=0
velFile.close()
swimFile.close()
if(fieldType=="nem"):
  nemFile.close()

# Animate
print( "Animating ..." )
if(fieldType=="vel"):
  name='%s/swimmerVelField_animation%s'%(dataPath,suffix)
elif(fieldType=="vor"):
  name='%s/swimmerVorField_animation%s'%(dataPath,suffix)
elif(fieldType=="nem"):
  name='%s/swimmerDirField_animation%s'%(dataPath,suffix)
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
  myCommand="rm frame*.png"
  call(myCommand,shell=True)
