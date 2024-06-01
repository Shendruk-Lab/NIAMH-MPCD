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
vel_vor = sys.argv[2] # string 

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

vel_vor.lower()
if not (vel_vor == 'vel' or vel_vor == 'vor'):
    print("has to be 'vel' or 'vor'")
    exit()

###########################################################
### Format and style
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
if (vel_vor == 'vel'):
	myMap=ed.deepsea
else:
    myMap=ed.bombpops
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
for i in range(21): 
    buff = file.readline()
buff = file.readline().split() #for the chunk colours
dStrength=float(buff[-1].replace('(','').replace(')',''))
buff = file.readline().split()
dChunks = int(buff[-1].replace('(','').replace(')',''))
for i in range(36):
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

# for chunk colours: red extensile, blue contractile, black neutral
lenChunk = int(monoN/dChunks)
maxN = lenChunk*dChunks
dipoles = []
for i in range(monoN):
    checker = (i//lenChunk)%2
    print(checker)
    if i<maxN:
        if checker == 0:
            dipoles.append(dStrength)
        else:
            dipoles.append(-dStrength)
    else:
        dipoles.append(0.0)
print("dipoles: ", dipoles)
    


MDMPCD=int(0.1/dtMD)                             # number of MD timesteps per MPCD timestep
dirSOut=int(outMD/MDMPCD)        				 # the frequency of mpcd output 

nemName=mpcdDataPath+'/velfield.dat'
sName=mpcdDataPath+'/directorfield.dat'
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
VEL = zeros(shape=(3,xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
V = zeros(shape=(2,xyzSize[0],xyzSize[1]),dtype=float)
S = zeros(shape=(xyzSize[0],xyzSize[1]),dtype=float)
X = zeros(shape=(xyzSize[0],xyzSize[1]),dtype=float)
Y = zeros(shape=(xyzSize[0],xyzSize[1]),dtype=float)

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
# print( "\tReanalysisScripts/polymer2Dplot-edited2.pyad frame %d"%(n+1) )
while nemFile:
	i=i+1
	# Read nematic field
	line = nemFile.readline()
	if( not line ):
		break
	else:
		if j>=start:
			t,Qx,Qy,Qz,Vx,Vy,Vz = line.split("\t",7)
			XYZ[0][int(Qx)][int(Qy)][int(Qz)] = float(Qx)
			XYZ[1][int(Qx)][int(Qy)][int(Qz)] = float(Qy)
			XYZ[2][int(Qx)][int(Qy)][int(Qz)] = float(Qz)
			VEL[0][int(Qx)][int(Qy)][int(Qz)] = float(Vx)
			VEL[1][int(Qx)][int(Qy)][int(Qz)] = float(Vy)
			VEL[2][int(Qx)][int(Qy)][int(Qz)] = float(Vz)
			# S[int(Qx)][int(Qy)][int(Qz)] = float(s)
			#S[int(Qx)][int(Qy)] = float(s)
	
	if i==xyzSize[0]*xyzSize[1]*xyzSize[2]:
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
		# Calculate speed/vorticity
		z=int(cm[2])
		for x in range(xyzSize[0]):
			for y in range(xyzSize[1]):
				X[x][y]=XYZ[0][x][y][z]
				Y[x][y]=XYZ[1][x][y][z]
				for d in range(2):
					V[d][x][y]=VEL[d][x][y][z]
		if(vel_vor=='vel'):
			for x in range(xyzSize[0]):
				for y in range(xyzSize[1]):
					S[x][y]=sqrt(VEL[0][x][y][z]**2+VEL[1][x][y][z]**2+VEL[2][x][y][z]**2)
		if(vel_vor=='vor'): #stolenish from swimmerTraj_fields2D.py
			dx=1.0
			dy=1.0
			# corners both x and y are forward or reverse derivatives
			#0,0
			vydx=(V[1][1][0]-V[1][0][0])/dx
			vxdy=(V[0][0][1]-V[0][0][0])/dy
			S[0][0]=vydx-vxdy
			#N,0
			vydx=(V[1][-1][0]-V[1][-2][0])/dx
			vxdy=(V[0][-1][0]-V[0][-1][1])/dy
			S[-1][0]=vydx-vxdy
			#N,N
			vydx=(V[1][-1][-1]-V[1][-2][-1])/dx
			vxdy=(V[0][-1][-1]-V[0][-1][-2])/dy
			S[-1][-1]=vydx-vxdy
			#0,N
			vydx=(V[1][1][-1]-V[1][0][-1])/dx
			vxdy=(V[0][0][-1]-V[0][0][-2])/dy
			S[0][-1]=vydx-vxdy
   
			# edges are forward or reverse with centred derivatives
			for x in range(1, xyzSize[0]-1):
				vydx=0.5*(V[1][x+1][0]-V[1][x-1][0])/dx
				vxdy=(V[0][x][-1]-V[0][x][0])
				vydx=0.5*(V[1][x+1][-1]-V[1][x-1][-1])/dx
				vxdy=(V[0][x][-1]-V[0][x][-2])/dy
				S[x][-1] =vydx-vxdy
			for y in range(1,xyzSize[1]-1):
				vydx=(V[1][1][y]-V[1][0][y])/dx
				vxdy=0.5*(V[0][0][y+1]-V[0][0][y-1])/dy
				S[0][y] =vydx-vxdy
				vydx=(V[1][-1][y]-V[1][-2][y])/dx
				vxdy=0.5*(V[0][-1][y+1]-V[0][-1][y-1])/dy
				S[-1][y] =vydx-vxdy
   
			#The bulk is all centred derivatives
			for x in range(1,xyzSize[0]-1):
				for y in range(1,xyzSize[1]-1):
					vydx=0.5*(V[1][x+1][y]-V[1][x-1][y])/dx
					vxdy=0.5*(V[0][x][y+1]-V[0][x][y-1])/dy
					S[x][y] =vydx-vxdy
		j=j+1
		if j>finish:
			break
		# Save the instantaneous or current velocity field frame
		# Save frame
		n=n+1
		print( "\t\tPlot frame %d"%(n) )
		# Setup figure object
		tight_layout()
		fig,ax = plt.subplots()
		plt.cla()
		if(projDim==2):
			quiver( X,Y,V[0],V[1] )
			CS3 = imshow(S.T,cmap=myMap,origin='lower')
			cb = fig.colorbar(CS3)
			if(vel_vor=='vel'):
				cb.ax.set_ylabel(r'$v$', fontsize = FS)
			else:
				cb.ax.set_ylabel(r'$\omega$', fontsize = FS)
			# circle = Circle((cm[j-1][0],cm[j-1][1]),radius = 0.5, facecolor = "none", edgecolor = "black",zorder =1)
			# if (j-1>0):
			# 	for k in range(j-1):
			# 		plot( [cm[k][0],cm[k+1][0]],[cm[k][1],cm[k+1][1]],color='green',linewidth=2,zorder=2 )
			for k in range(monoN-1):
				if(fabs(POLY[0][k]-POLY[0][k+1])<0.5*xyzSize[0] and fabs(POLY[1][k]-POLY[1][k+1])<0.5*xyzSize[1]):
					plot( [POLY[0][k],POLY[0][k+1]],[POLY[1][k],POLY[1][k+1]],color='black',linewidth=2,zorder=2 )
			for k in range(monoN):
				if(dipoles[k]>0.000001):
					if(vel_vor=='vel'):
						plot([POLY[0][k]],[POLY[1][k]], 'o', color='tab:orange', markersize=5)
					else:
						plot([POLY[0][k]],[POLY[1][k]], 'o', color='r', markersize=5)
				elif(dipoles[k]<-0.000001):
					if(vel_vor=='vel'):
						plot([POLY[0][k]],[POLY[1][k]], 'o', color='m', markersize=5)
					else:
						plot([POLY[0][k]],[POLY[1][k]], 'o', color='b', markersize=5)
				else:
					plot([POLY[0][k]],[POLY[1][k]], 'ko', markersize=5)
			fig.canvas.draw()
		"""
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
						yline=[ XYZ[1][x][y][z]-c*VEL[1][x][y][z],XYZ[1][x][y][z]+c*VEL[1][x][y][z] ]
						zline=[ XYZ[2][x][y][z]-c*VEL[2][x][y][z],XYZ[2][x][y][z]+c*VEL[2][x][y][z] ]
						plot( yline,zline,color=myMap(S[x][y][z],1),linewidth=myLW )
			for k in range(monoN-1):
				if(fabs(POLY[1][k]-POLY[1][k+1])<0.5*xyzSize[1] and fabs(POLY[2][k]-POLY[2][k+1])<0.5*xyzSize[2]):
					if(int(dipoles[k]+dipoles[k+1])>0): #extensile
						print("extensile")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='r',linewidth=2 )
					elif(int(dipoles[k]+dipoles[k+1])<0): #contractile
						print("contractile")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='b',linewidth=2 )
					else:
						print("link")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='k',linewidth=2 )
		elif(projDim==1):
			y=int(cm[1])
			for x in range(xyzSize[0]):
				for z in range(xyzSize[2]):
					if( S[x][y][z]>0.05 ):
						xline=[ XYZ[0][x][y][z]-c*VEL[0][x][y][z],XYZ[0][x][y][z]+c*VEL[0][x][y][z] ]
						zline=[ XYZ[2][x][y][z]-c*VEL[2][x][y][z],XYZ[2][x][y][z]+c*VEL[2][x][y][z] ]
						plot( xline,zline,color=myMap(S[x][y][z],1),linewidth=myLW )
			for k in range(monoN-1):
				if(fabs(POLY[0][k]-POLY[0][k+1])<0.5*xyzSize[0] and fabs(POLY[2][k]-POLY[2][k+1])<0.5*xyzSize[2]):
					if(int(dipoles[k]+dipoles[k+1])>0): #extensile
						print("extensile")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='r',linewidth=2 )
					elif(int(dipoles[k]+dipoles[k+1])<0): #contractile
						print("contractile")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='b',linewidth=2 )
					else:
						print("link")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='k',linewidth=2 )
		elif(projDim==2):
			print("dimensions 2")
			z=int(cm[2])
			for x in range(xyzSize[0]):
				for y in range(xyzSize[1]):
					# if( S[x][y][z]>0.05 ):
					if( S[x][y]>0.05 ):
						xline=[ XYZ[0][x][y][z]-c*VEL[0][x][y][z],XYZ[0][x][y][z]+c*VEL[0][x][y][z] ]
						yline=[ XYZ[1][x][y][z]-c*VEL[1][x][y][z],XYZ[1][x][y][z]+c*VEL[1][x][y][z] ]
						# plot( xline,yline,color=myMap(S[x][y][z],1),linewidth=myLW )
						plot( xline,yline,color="white",linewidth=myLW )
			imshow(S.T,cmap=myMap,origin='lower',vmin=0,vmax=1)
			# circle = Circle((cm[j-1][0],cm[j-1][1]),radius = 0.5, facecolor = "none", edgecolor = "black",zorder =1)
			# if (j-1>0):
			# 	for k in range(j-1):
			# 		plot( [cm[k][0],cm[k+1][0]],[cm[k][1],cm[k+1][1]],color='green',linewidth=2,zorder=2 )
			for k in range(monoN-1):
				if(fabs(POLY[0][k]-POLY[0][k+1])<0.5*xyzSize[0] and fabs(POLY[1][k]-POLY[1][k+1])<0.5*xyzSize[1]):
					if(int(dipoles[k]+dipoles[k+1])>0): #extensile
						print("extensile")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='r',linewidth=2, zorder=2 )
					elif(int(dipoles[k]+dipoles[k+1])<0): #contractile
						print("contractile")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='b',linewidth=2, zorder=2 )
					else:
						print("link")
						plot( [POLY[1][k],POLY[1][k+1]],[POLY[2][k],POLY[2][k+1]],color='black',linewidth=2, zorder=2 )

		"""
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
		if(vel_vor=='vel'):
			name='velframe%04d.png'%(n)
		else:
			name='vorframe%04d.png'%(n)
		savefig( name, format='png', dpi=600, bbox_inches='tight')
		i=0

nemFile.close()
polyFile.close()

###########################################################
### Animate data
###########################################################
# print( "\tAnimating ..." )
if vel_vor=='vel':
	name='velocity_animation%s'%suffix
else:
    name='vorticity_animation%s'%suffix
#name='2DPolymerProj_animation%s'%suffix
myCommand="rm %s"%name
call(myCommand,shell=True)
if(vel_vor=='vel'):
    myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i velframe%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
else:
    myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i vorframe%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
	myCommand="rm frame*.png"
	call(myCommand,shell=True)
