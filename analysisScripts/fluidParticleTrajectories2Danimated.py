"""
    Fluid particle and mobile colloid trajectory rendering script.

    Originally from Tyler N. Shendruk
    Modified by Timofey Kozhukhov
"""

from pylab import *
from subprocess import call
import sys
import os
import json
import argparse
from matplotlib import patches

###########################################################
### Set up argsparse
###########################################################
parser = argparse.ArgumentParser(description='Particle trajectory rendering script.')
parser.add_argument('-i',"--input", help="Path to input .json file")
parser.add_argument('-p',"--particle", nargs='+', help='List of paths to the fluid particle data (should be "detailedSP0.dat")')
parser.add_argument('-c',"--colloid", nargs='+', help='List of paths to the colloid particle data (should be "solidtraj0.dat")')
parser.add_argument('-s',"--start", type=int, help="Starting timestep for averaging",default=0)
parser.add_argument('-f',"--finish", type=int, help="Finishing timestep for averaging",default=999999999)
parser.add_argument("-k", "--keepFrames", type=int,help="0=don't keep (delete) frames; 1=keep frames",default=0)
args = parser.parse_args()

###########################################################
### Read arguments
###########################################################
print("Arguments:")
for arg, value in vars(args).items():
    print(f"\t{arg}: {value}")
inputName = args.input
particleName = args.particle
colloidName = args.colloid
start = args.start
finish = args.finish
keepFrames = args.keepFrames

###########################################################
### Format and style
###########################################################
# Use our custom style and colours
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed
# Colour map to use
myMap=ed.plasma
# Adjust line width
myLW=3.0
# adjust line length
c = 1.0
#Animation stuff
bitrate=5000
framerate=12            #Number of frames per second in the output video
codec='libx264'         #Other options include mpeg4
suffix='.mp4'

###########################################################
### Initialize
###########################################################
_x=0
_y=1
_z=2

###########################################################
### Read input json
###########################################################
if not os.path.isfile(inputName):
    print("%s not found."%inputName)
    exit()
with open(inputName, 'r') as f:
  input = json.load(f)
xyzSize=array([30,30,1])
dim=2
if "domain" in input:
    xyzSize[_x]=input['domain'][_x]
    dim=1
    if(len(input['domain'])>1):
        xyzSize[_y]=input['domain'][_y]
        dim=2
    else:
        xyzSize[_y]=1
    if(len(input['domain'])>2):
        xyzSize[_z]=input['domain'][_z]
        dim=3
    else:
        xyzSize[_z]=1
LC=0
if "lc" in input:
    LC=input['lc']
popFluid=[]
if "species" in input:
    for s in input['species']:
        if 'pop' in s:
            popFluid.append(s['pop'])
        else:
            print("Unfortunately this script does not support species without 'pop' in input json.")
            exit()
radius=[]
aspectRatio=[]
if "BC" in input:
    for c in input['BC']:
        if 'R' in c:
            radius.append(c['R'])
        else:
            print("Unfortunately this script does not support BCs without 'R' in input json.")
            exit()
        if 'aInv' in c:
            aspectRatio.append(c['aInv'][_x]/c['aInv'][_y])
        else:
            aspectRatio.append(1.0)
solidTrajOut=0
trajOut=0
if "solidTrajOut" in input:
    solidTrajOut=input['solidTrajOut']
if "trajOut" in input:
    trajOut=input['trajOut']

###########################################################
### Check input and arguments match
###########################################################
if(len(particleName) != len(popFluid)):
    print("Number of particle data files does not match number of species in input json.")
    exit()
if(len(colloidName) != len(radius)):
    print("Number of colloid data files does not match number of BCs in input json.")
    exit()
if(len(colloidName) != len(aspectRatio)):
    print("Number of colloid data files does not match number of aspect ratios in input json.")
    exit()
numFluidParticles = len(particleName)
numColloids = len(colloidName)

###########################################################
### Setup the animation
###########################################################
# Figure
fig = plt.figure(1)
ax = fig.add_subplot(111, aspect='equal')

###########################################################
### Read the data for animation
###########################################################
print( 'Reading data ...' )
print( '\tReading colloid data ...' )
cPos=[]
cOri=[]
for c in colloidName:
    print( '\t\tRead %s.'%c )
    if not os.path.isfile(c):
        print("%s not found."%c)
        exit()
    infile = open(c,"r")
    # Toss header
    for i in range(13):
        line = infile.readline()
    i=0
    cPos.append([])
    cOri.append([])
    while infile:
        line = infile.readline()
        if( not line ):
            break
        else:
            if i>=start:
                t,Qx,Qy,Qz,_,_,_,Ox,Oy,Oz,_,_,_ = line.split("\t",13)
                cPos[-1].append([float(Qx),float(Qy)])
                cOri[-1].append(float(Oz))
        i=i+1
    infile.close()
print( '\tReading fluid particle data ...' )
fPos=[]
fOri=[]
for c in range(numFluidParticles):
    p=particleName[c]
    print( '\t\tRead %s.'%p )
    if not os.path.isfile(p):
        print("%s not found."%p)
        exit()
    infile = open(p,"r")
    # Toss header
    for i in range(14):
        line = infile.readline()
    i=0
    fPos.append([])
    fOri.append([])
    while infile:
        fPos[c].append([[],[]])
        fOri[c].append([])
        for j in range(popFluid[c]):
            line = infile.readline()
            if( not line ):
                break
            else:
                if i>=start:
                    t,Qx,Qy,Qz,_,_,_,_,Ux,Uy,Uz = line.split("\t",11)
                    fPos[c][i][_x].append(float(Qx))
                    fPos[c][i][_y].append(float(Qy))
                    fOri[c][i].append([float(Ux),float(Uy)])
        line = infile.readline()
        if( not line ):
            break
        i=i+1
    infile.close()
time=len(cPos[0])
for i in range(len(fPos)):
    time=min(time,len(fPos[i]))
for i in range(len(cPos)):
    time=min(time,len(cPos[i]))

print(inputName, particleName, colloidName, start, finish, keepFrames)
print(xyzSize,LC,popFluid, radius, aspectRatio)
print(solidTrajOut, trajOut)

for t in range(time):
    fig = plt.figure(1)
    plt.cla()
    for c in range(numFluidParticles):
        if(not LC):
            scatter(fPos[c][t][_x],fPos[c][t][_y])
        else:
            for i in range(len(fOri[c][t])):
                if(i==0):
                    plot([fPos[c][t][_x][i]-0.5*fOri[c][t][i][_x],fPos[c][t][_x][i]+0.5*fOri[c][t][i][_x]],[fPos[c][t][_y][i]-0.5*fOri[c][t][i][_y],fPos[c][t][_y][i]+0.5*fOri[c][t][i][_y]],linewidth=myLW)
                    colour=ax.get_lines()[-1].get_c()
                else:
                    plot([fPos[c][t][_x][i]-0.5*fOri[c][t][i][_x],fPos[c][t][_x][i]+0.5*fOri[c][t][i][_x]],[fPos[c][t][_y][i]-0.5*fOri[c][t][i][_y],fPos[c][t][_y][i]+0.5*fOri[c][t][i][_y]],linewidth=myLW,color=colour)

    for c in range(numColloids):
        for i in [-1,0,1]:
            for j in [-1,0,1]:
                ell = patches.Ellipse((i*xyzSize[_x]+cPos[c][t][_x],j*xyzSize[_y]+cPos[c][t][_y]),2*radius[c]*aspectRatio[c],2*radius[c],angle=cOri[c][t]*180/np.pi,linewidth=2, fill=True, zorder=2, color='k')
                ax.add_patch(ell)
    xlabel(r'$x$')
    ylabel(r'$y$')
    plt.axis(xmax=xyzSize[_x],xmin=0,ymax=xyzSize[_y],ymin=0)
    name='frame%04d.png'%(t)
    # uncomment below for snapshots
    plt.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    savefig( name,bbox_inches='tight',pad_inches=0 )

###########################################################
### Visualize data
###########################################################
print( "Animating ..." )
name='2Dorientation_animation%s'%suffix
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)
if not keepFrames:
    myCommand="rm frame*.png"
    call(myCommand,shell=True)
