###########################################################
### Animates 2D fields in swimmers' reference frame
###########################################################

###########################################################
### Imports
###########################################################
from pylab import *
from subprocess import call

###########################################################
### Read arguments
###########################################################
FS =25
TLS = 20		# Tick label size
xyzSize=zeros( 3,dtype=int )
print( "Arguments:" )
for arg in sys.argv:
    print( "\t" + arg )
dataPath = sys.argv[1]		# Name of the data
numSw = int(sys.argv[2])	# Number of swimmers --- reference frame shifted about FIRST swimmer
xyzSize[0] = int(sys.argv[3])	# System size
xyzSize[1] = int(sys.argv[4])	# System size
xyzSize[2] = int(sys.argv[5])	# System size
dipole = float(sys.argv[6])	# Dipole size
start = int(sys.argv[7])		# Average after this number
finish = int(sys.argv[8])		# Average before this number
qx = int(sys.argv[9])		# Only show every qx arrow in x
qy = int(sys.argv[10])		# Only show every qy arrow in y
avdim = sys.argv[11]			# Dimension to average over
rotAx=int(sys.argv[12])		#0=just shift the swimmer to the centre; 1=it aligned in x-direction too
myAspect=sys.argv[13]		#'auto' - reshapes into square graph or 'equal' keeps whatever aspect ratio the true values
makeAni=int(sys.argv[14])	#0=just make the average (no animation); 1=make animation too
keepFrames=int(sys.argv[15])	#0=don't keep (delete) frames; 1=keep frames

###########################################################
### Set arguments
###########################################################
tailFreq=1.0/3.5
tailWaveLength=(4.0/10.0)*2.0*(1.0+dipole)
hidgeonLength=1.0
tailRad=1.0

###########################################################
### Style/formating stuff
###########################################################
# Style
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed
# Colours
swimmerMap=ed.viridis
myMap = ed.truncate_colormap(ed.viridis, minval=0.1, maxval=0.7)
#Animation
bitrate=5000
framerate=20		#Number of frames per second in the output video
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

###########################################################
### Initialize
###########################################################
velFieldName="%s/flowfield.dat"%dataPath
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
### Setup the animation
fig1, axes = plt.subplots(nrows=1, ncols=1)
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
# Data
XYZ = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
VEL = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
MAG = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
# Shifted/rotated frames
rotSize = int(sqrt(xyzSize[d1]**2+xyzSize[d2]**2))+1
shiftMEAN = zeros(shape=(3,rotSize,rotSize),dtype=float)
shiftMAG = zeros(shape=(rotSize,rotSize),dtype=float)
rotMEAN = zeros(shape=(3,rotSize,rotSize),dtype=float)
rotMAG = zeros(shape=(rotSize,rotSize),dtype=float)
avRotMEAN = zeros(shape=(3,rotSize,rotSize),dtype=float)
avRotMAG = zeros(shape=(rotSize,rotSize),dtype=float)
cnt=zeros(shape=(rotSize,rotSize),dtype=float)
# Make Mesh
XY = zeros(shape=(2,rotSize,rotSize),dtype=float)
for x in range(rotSize):
    for y in range(rotSize):
        XY[0][x][y]=x
        XY[1][x][y]=y
###########################################################
### Read the data for min/max
###########################################################
print( 'Read file for min/max ...' )
velFile = open(velFieldName,"r")
swimFile = open(swimmerName,"r")
minV=99999999999999.0
maxV=0.0

###########################################################
### Read the data for animation
###########################################################
print( 'Read files for animation ...' )
print( '\tToss headers ...' )
for i in range(13):
  line = velFile.readline()
for i in range(13):
  line = swimFile.readline()
print( '\tRead data ...' )
i=0
j=0
n=-1
while velFile:
  i=i+1
  line = velFile.readline()
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
    ###########################################################
    ### Average the velocity
    ###########################################################
    j=j+1
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
        for y in range(xyzSize[1]):
          for z in range(xyzSize[2]):
            for d in range(3):
              MEAN[d][y][z]/=xyzSize[0]
      elif avdim=='y':
        for x in range(xyzSize[0]):
          for y in range(xyzSize[1]):
            for z in range(xyzSize[2]):
              MEAN[0][x][z]=MEAN[0][x][z]+VEL[0][x][y][z]
              MEAN[1][x][z]=MEAN[1][x][z]+VEL[1][x][y][z]
              MEAN[2][x][z]=MEAN[2][x][z]+VEL[2][x][y][z]
        for x in range(xyzSize[0]):
          for z in range(xyzSize[2]):
            for d in range(3):
                MEAN[d][x][z]/=xyzSize[1]
      elif avdim=='z':
        for x in range(xyzSize[0]):
          for y in range(xyzSize[1]):
            for z in range(xyzSize[2]):
              MEAN[0][x][y]=MEAN[0][x][y]+VEL[0][x][y][z]
              MEAN[1][x][y]=MEAN[1][x][y]+VEL[1][x][y][z]
              MEAN[2][x][y]=MEAN[2][x][y]+VEL[2][x][y][z]
        for x in range(xyzSize[0]):
          for y in range(xyzSize[1]):
            for d in range(3):
              MEAN[d][x][y]/=xyzSize[2]
        for x in range(xyzSize[d1]):
          for y in range(xyzSize[d2]):
            MAG[x][y]=sqrt( MEAN[0][x][y]**2+MEAN[1][x][y]**2+MEAN[2][x][y]**2 )
        #Save the instantaneous or current velocity field frame
        for x in range(xyzSize[d1]):
            for y in range(xyzSize[d2]):
                if MAG[x][y]>maxV:
                    maxV=MAG[x][y]
                elif MAG[x][y]<minV:
                    minV=MAG[x][y]
    ###########################################################
    ### Read the swimmers' positions
    ###########################################################
    H=[]  #Head
    B=[]  #Body
    T=[]  #Tail
    for ns in range(numSw):
        line = swimFile.readline()
        if( not line ):
            break
        else:
            t,Hx,Hy,Hz,Hvx,Hvy,Hvz,Tx,Ty,Tz,Tvx,Tvy,Tvz,runTumble = line.split( )
            H.append( [float(Hx),float(Hy),float(Hz)] )
            B.append( [float(Tx),float(Ty),float(Tz)] )
            T.append( [(1.0+dipole)*B[-1][0]-dipole*H[-1][0],(1.0+dipole)*B[-1][1]-dipole*H[-1][1],(1.0+dipole)*B[-1][2]-dipole*H[-1][2]] )
    #Toss the line break
    line = swimFile.readline()
    if( not line ):
        break
    ###########################################################
    ### Shift everything to be about the CM of the "first" swimmer
    ###########################################################
    shift=[]
    for d in range(3):
        cm=0.5*(H[0][d]+T[0][d])
        shift.append( 0.5*float(rotSize)-cm )
    if avdim=='x':
        shift[0]=0.0
    elif avdim=='y':
        shift[1]=0.0
    elif avdim=='z':
        shift[2]=0.0
    # Shift the swimmers
    for ns in range(numSw):
        for d in range(3):
            H[ns][d]+=shift[d]
            B[ns][d]+=shift[d]
            T[ns][d]+=shift[d]
    # Shift the flow field
    for x in range(xyzSize[d1]):
        xshift=x+shift[d1]
        if( xshift>0.5*(rotSize+xyzSize[d1]) ):
            xshift-=xyzSize[d1]
        elif( xshift<0.5*(rotSize-xyzSize[d1]) ):
            xshift+=xyzSize[d1]
        xshift=int(xshift)
        for y in range(xyzSize[d2]):
            yshift=y+shift[d2]
            if( yshift>0.5*(rotSize+xyzSize[d2]) ):
                yshift-=xyzSize[d2]
            elif( yshift<0.5*(rotSize-xyzSize[d2]) ):
                yshift+=xyzSize[d2]
            yshift=int(yshift)
            for d in range(3):
                shiftMEAN[d][xshift][yshift]=MEAN[d][x][y]
            shiftMAG[xshift][yshift]=MAG[x][y]
    ###########################################################
    ### Rotate everything so that the first swimmer is along x-dir
    ###########################################################
    if rotAx:
        theta=-arctan2( (H[0][d2]-T[0][d2]), (H[0][d1]-T[0][d1]) )
        R=array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])
        cm=[ 0.5*(H[0][d1]+T[0][d1]),0.5*(H[0][d2]+T[0][d2]) ]
        # Rotate the swimmers
        for ns in range(numSw):
            x=(H[ns][d1]-cm[0])*R[0][0]+(H[ns][d2]-cm[1])*R[0][1]
            y=(H[ns][d1]-cm[0])*R[1][0]+(H[ns][d2]-cm[1])*R[1][1]
            H[ns][d1]=x+cm[0]
            H[ns][d2]=y+cm[1]
            x=(B[ns][d1]-cm[0])*R[0][0]+(B[ns][d2]-cm[1])*R[0][1]
            y=(B[ns][d1]-cm[0])*R[1][0]+(B[ns][d2]-cm[1])*R[1][1]
            B[ns][d1]=x+cm[0]
            B[ns][d2]=y+cm[1]
            x=(T[ns][d1]-cm[0])*R[0][0]+(T[ns][d2]-cm[1])*R[0][1]
            y=(T[ns][d1]-cm[0])*R[1][0]+(T[ns][d2]-cm[1])*R[1][1]
            T[ns][d1]=x+cm[0]
            T[ns][d2]=y+cm[1]
        # Rotate the flow field
        for x in range(rotSize):
            for y in range(rotSize):
                xrot=(x-cm[0])*R[0][0]+(y-cm[1])*R[0][1]
                yrot=(x-cm[0])*R[1][0]+(y-cm[1])*R[1][1]
                xrot=int(xrot+cm[0])
                yrot=int(yrot+cm[1])
                if(xrot>=rotSize):
                    xrot-=rotSize
                elif(xshift<0):
                    xrot+=rotSize
                if(yrot>=rotSize):
                    yrot-=rotSize
                elif(yrot<0):
                    yrot+=rotSize
                # Rotate the actual velocity vectors
                if(shiftMAG[x][y]>0.0):
                    rotMEAN[d1][xrot][yrot]=shiftMEAN[d1][x][y]*R[0][0]+shiftMEAN[d2][x][y]*R[0][1]
                    rotMEAN[d2][xrot][yrot]=shiftMEAN[d1][x][y]*R[1][0]+shiftMEAN[d2][x][y]*R[1][1]
                    rotMAG[xrot][yrot]=shiftMAG[x][y]
    else:
        for x in range(rotSize):
            for y in range(rotSize):
                rotMEAN[d1][x][y]=shiftMEAN[d1][x][y]
                rotMEAN[d2][x][y]=shiftMEAN[d2][x][y]
                rotMAG[x][y]=shiftMAG[x][y]
    ###########################################################
    ### Average the shifted and rotate field
    ###########################################################
    for x in range(rotSize):
        for y in range(rotSize):
            if(rotMAG[x][y]>0.0):
                cnt[x][y]+=1.0
                avRotMAG[x][y]+=rotMAG[x][y]
                for d in range(3):
                    avRotMEAN[d][x][y]+=rotMEAN[d][x][y]
    ###########################################################
    ### Plot the frame
    ###########################################################
    # Save frame
    n=n+1
    if(makeAni):
      #Setup the density image
      plt.subplot(1,1,1)
      #Setup the velocity image
      quiv = quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], rotMEAN[d1][::qx, ::qy], rotMEAN[d2][::qx, ::qy] )
      velImage = imshow(rotMAG.T,cmap=myMap,origin='lower',aspect=myAspect,vmin=minV,vmax=maxV)
      velCB = colorbar()
      velCB.ax.set_ylabel(r'Velocity, $\left|\vec{v}\right|$', fontsize = FS)
      # Plot the swimmers
      for ns in range(numSw):
        if( fabs(H[ns][d1]-B[ns][d1])<0.5*xyzSize[d1] and fabs(H[ns][d2]-B[ns][d2])<0.5*xyzSize[d2] ):
          plot( [H[ns][d1],B[ns][d1]],[H[ns][d2],B[ns][d2]],'-',color='k',linewidth=2 )
        if( fabs(T[ns][d1]-B[ns][d1])<0.5*xyzSize[d1] and fabs(T[ns][d2]-B[ns][d2])<0.5*xyzSize[d2] ):
          plot( [T[ns][d1],B[ns][d1]],[T[ns][d2],B[ns][d2]],'-',color=swimmerMap(B[ns][dim]/xyzF[dim]),linewidth=2 )
        waveLength=(4.0/5.0)*2.0*(1.0+dipole)
        theta=arctan2(H[ns][d2]-B[ns][d2],H[ns][d1]-B[ns][d1])
        X=linspace(0,2.0*(1.0+dipole),20)
        # Y=sin(2.0*pi*X/waveLength)*sin(0.5*pi + float(n)*tailFreq)
        Y=tailRad*(1.0-exp(-pow(X/hidgeonLength,2)))*sin(2.0*pi*X/tailWaveLength + float(n)*tailFreq)
        tailX=cos(theta)*X - sin(theta)*Y
        tailY=sin(theta)*X + cos(theta)*Y
        plot( B[ns][d1]-tailX,B[ns][d2]-tailY,'-',color=swimmerMap(B[ns][dim]/xyzF[dim]),linewidth=2 )
        plot( H[ns][d1],H[ns][d2],'o',color='k',fillstyle='full' )
        plot( B[ns][d1],B[ns][d2],'h',color=swimmerMap(B[ns][dim]/xyzF[dim]),fillstyle='full' )
      xlabel(r'$%s$'%labX, fontsize = FS)
      ylabel(r'$%s$'%labY, fontsize = FS)
      #plt.axis(xmax=xyzSize[0], xmin=0, ymax=xyzSize[1], ymin=0)
      plt.axis(xmax=rotSize, xmin=0, ymax=rotSize, ymin=0)
      name='frame%04d.png'%(n)
      savefig( name )
      plt.clf()
    #Zero matrix
    VEL= zeros( (3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float )
    MEAN = zeros(shape=(3,xyzSize[d1],xyzSize[d2]),dtype=float)
    shiftMEAN = zeros(shape=(3,rotSize,rotSize),dtype=float)
    shiftMAG = zeros(shape=(rotSize,rotSize),dtype=float)
    rotMEAN = zeros(shape=(3,rotSize,rotSize),dtype=float)
    rotMAG = zeros(shape=(rotSize,rotSize),dtype=float)
    i=0
velFile.close()
swimFile.close()

###########################################################
### Animate
###########################################################
if(makeAni):
    print( "\tAnimating ..." )
    plt.close("all")
    if rotAx:
        name='%s/swimmerVelField_refFrameRotated_animation%s'%(dataPath,suffix)
    else:
        name='%s/swimmerVelField_RefFrameCentred_animation%s'%(dataPath,suffix)
    myCommand="rm %s"%name
    call(myCommand,shell=True)
    myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
    call(myCommand,shell=True)
    if not keepFrames:
        myCommand="rm *.png"
        call(myCommand,shell=True)

###########################################################
### Plot the average
###########################################################
minV=99999999999999.0
maxV=0.0
for x in range(rotSize):
    for y in range(rotSize):
        if(cnt[x][y]>0.0):
            avRotMAG[x][y]/=cnt[x][y]
            if(avRotMAG[x][y]<minV):
                minV=avRotMAG[x][y]
            elif(avRotMAG[x][y]>maxV):
                maxV=avRotMAG[x][y]
            for d in range(3):
                avRotMEAN[d][x][y]/=cnt[x][y]
#Setup the density image
plt.subplot(1,1,1)
plt.cla()
#Setup the velocity image
quiv = quiver( XY[0][::qx, ::qy], XY[1][::qx, ::qy], avRotMEAN[d1][::qx, ::qy], avRotMEAN[d2][::qx, ::qy] )
velImage = imshow(avRotMAG.T,cmap=myMap,origin='lower',aspect=myAspect,vmin=minV,vmax=maxV)
if(not makeAni):
    velCB = colorbar()
    velCB.ax.set_ylabel(r'Velocity, $\left|\vec{v}\right|$', fontsize = FS)
xlabel(r'$%s$'%labX, fontsize = FS)
ylabel(r'$%s$'%labY, fontsize = FS)
#plt.axis(xmax=xyzSize[0], xmin=0, ymax=xyzSize[1], ymin=0)
plt.axis(xmax=rotSize, xmin=0, ymax=rotSize, ymin=0)
if rotAx:
    name='%s/swimmerVelField_refFrameRotated_av.pdf'%(dataPath)
else:
    name='%s/swimmerVelField_RefFrameCentred_av.pdf'%(dataPath)
savefig( name )
show()
