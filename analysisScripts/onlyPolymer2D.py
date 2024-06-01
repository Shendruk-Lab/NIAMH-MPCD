from matplotlib import pyplot as plt
from pylab import *
from numpy import ma
from subprocess import call
import os
import sys

def func(dir):
    rcParams["figure.autolayout"] = True
    ###########################################################
    ### Plots 2D averaging over user defined direction
    ###########################################################

    ###########################################################
    ### Read arguments
    ###########################################################

    xyzSize=zeros( 3,dtype=int )
    #print( "Arguments:" )
    #for arg in sys.argv:
    #	print( "\t" + arg )
    mpcdDataPath = dir#sys.argv[1]		# Path to directory containing MPCD data (input.json, md.inp, directorfield.dat, folder for vmd file)
    ###########################################################
    ### Initialize
    ###########################################################
    projection = 'z'		# The direction to project over
    c = 0.5		#Length of director lines approx 0.5
    keepFrames=0	#0=don't keep (delete) frames; 1=keep frames
    start = 0		# Show after this number (in MPCD)
    finish = 999999999		# Finish after this number

    ############################################################
    qx = 0.5
    qy = 0.5

    ###########################################################
    ### Format and style
    ###########################################################
    # Use our custom style
    plt.style.use('shendrukGroupStyle')
    # Use our custom colours
    import shendrukGroupFormat as ed
    # Colour map to use
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
    #for l in range(3500):
    #    for k in range(monoN+2):
    #        line=polyFile.readline()
    #        if (not line):
    #            print("run out already")
    #            break
    n=-1
    while polyFile:
        line = polyFile.readline()
        if( not line ):
            break
        #if n>100:
        #    break
        line = polyFile.readline()
        # print(line)
        for k in range(monoN):
            line = polyFile.readline()
            # print line
            t,Qx,Qy,Qz = line.split()
            POLY[0][k]=float(Qx)
            POLY[1][k]=float(Qy)
            POLY[2][k]=float(Qz)
            
        #for k in range(monoN-1):
        #    if(fabs(POLY[0][k]-POLY[0][k+1])<0.5*xyzSize[0] and fabs(POLY[1][k]-POLY[1][k+1])<0.5*xyzSize[1]):
        #        plot( [POLY[0][k],POLY[0][k+1]],[POLY[1][k],POLY[1][k+1]],color='black',linewidth=2,zorder=2 )
        n=n+1
        print( "\t\tPlot frame %d"%(n) )
        # Setup figure object
        tight_layout()
        fig,ax = plt.subplots()
        plt.cla()
        CS3 = imshow(zeros(shape=(xySize[0],xySize[1]),dtype=float),cmap=myMap,origin='lower')
        #cb = fig.colorbar(CS3)
        for k in range(monoN-1):
            if(fabs(POLY[0][k]-POLY[0][k+1])<0.5*xyzSize[0] and fabs(POLY[1][k]-POLY[1][k+1])<0.5*xyzSize[1]):
                plot( [POLY[0][k],POLY[0][k+1]],[POLY[1][k],POLY[1][k+1]],color='black',linewidth=2,zorder=2 )
        for k in range(monoN):
            if(dipoles[k]>0.000001):
                plot([POLY[0][k]],[POLY[1][k]], 'o', color='r', markersize=2.5)
            elif(dipoles[k]<-0.000001):
                plot([POLY[0][k]],[POLY[1][k]], 'o', color='b', markersize=2.5)
            else:
                plot([POLY[0][k]],[POLY[1][k]], 'ko', markersize=2.5)
        #for k in range((monoN+2)*3):
        #    line = polyFile.readline()
        
        axis('off')
        grid(False)
        subplots_adjust(left=0, bottom=0, right=1, top=1)
        name='frame%04d.png'%(n)
        savefig( name, format='png', dpi=600, bbox_inches='tight')
        i=0

    polyFile.close()
    print("I have closed the file")
    ###########################################################
    ### Animate data
    ###########################################################
    # print( "\tAnimating ..." )
    name='2DPolymerProj_animation%s'%suffix
    myCommand="rm %s"%name
    call(myCommand,shell=True)
    myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
    call(myCommand,shell=True)
    if not keepFrames:
        myCommand="rm frame*.png"
        call(myCommand,shell=True)
        
func(os.getcwd())