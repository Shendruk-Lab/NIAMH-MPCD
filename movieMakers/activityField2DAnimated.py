from sunau import AUDIO_FILE_ENCODING_LINEAR_24
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
dataName = sys.argv[1]			# Name of the data (should be coarsegrain.dat)
xyzSize[0] = int(sys.argv[2])	# System size
xyzSize[1] = int(sys.argv[3])	# System size
xyzSize[2] = int(sys.argv[4])	# System size
start = int(sys.argv[5])		# Average after this number
finish = int(sys.argv[6])		# Average before this number
sliceDim = sys.argv[7]				# Dimension to average over
sliceIndex = int(sys.argv[8])	# Plane to be sliced
myAspect=sys.argv[9]			#'auto' - reshapes into square graph or 'equal' keeps whatever aspect ratio the true values
keepFrames=int(sys.argv[10])	#0=don't keep (delete) frames; 1=keep frames
activityType = int(sys.argv[11]) # -1=all inc density, 0=actSum, 1=actAv, 2=actSigmoid
particleActivity = float(sys.argv[12]) # activity of individual particles
# for actSig
actSigAv = 20
actSigWidth = 1
if activityType == 2:
    actSigAv = int(sys.argv[13]) # average particle per cell count
    actSigWidth = float(sys.argv[14]) # sigmoid width tuning parameter



###########################################################
### Format and style
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed

# Colour map to use 
myMap=ed.viridis 

#Animation stuff
bitrate=5000
framerate=12		#Number of frames per second in the output video
    # Note that the ideal for this is _very_ variable so play around with it
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

###########################################################
### Initialize
###########################################################
if sliceDim=='x':
    print("Note, for now we assume sliceDim=z! Edit the script to fix this as necessary")
    dim=0
    d1=1
    d2=2
elif sliceDim=='y':
    print("Note, for now we assume sliceDim=z! Edit the script to fix this as necessary")
    dim=1
    d1=0
    d2=2
elif sliceDim=='z':
    dim=2
    d1=0
    d2=1
else:
    print( "sliceDim must be 'x', 'y' or 'z' - not %s"%sliceDim )
    exit()
if sliceIndex>=xyzSize[dim]:
    print( "sliceIndex must be within the system: sliceIndex=%d but %s system size=%d"%(sliceIndex,sliceDim,xyzSize[dim]) )
    exit()

# Data
TOTPOP = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
sliceTOTPOP = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)

POP0, slicePop0, POP1, slicePOP1, slicePOP0prop, sliceRGBA = None, None, None, None, None, None

# Figure
fig1 = plt.figure(1)
###########################################################
### Read the data for animation
###########################################################

### Setup the animation
# Make labels
if sliceDim=='x':
    labX='y'
    labY='z'
elif sliceDim=='y':
    labX='x'
    labY='z'
elif sliceDim=='z':
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
i=0 # line counter (within timestep)
currTStep=0
n=-1
tStepMaxPop = 0
while infile:
    i=i+1
    line = infile.readline()
    if( not line ):
        break
    else:
        if currTStep>=start:
            currLine = line.split("\t")
            t,Qx,Qy,Qz,Vcmx,Vcmy,Vcmz,pop = [currLine[i] for i in range(8)] # get core data
            TOTPOP[int(Qx)][int(Qy)][int(Qz)] = float(pop)

    if i==xyzSize[0]*xyzSize[1]*xyzSize[2]: # end of a timestep
        currTStep=currTStep+1
        if currTStep>finish:
            break
        if currTStep<start or currTStep>finish:
      #print( 'Toss %d'%currTStep )
            aaa=0
        else:
            print( 'Work %d'%currTStep )

            # compute the information for the slice we want to render
            ##NOTE: for now, assuming sliceDim=z and z=0 for simplicity
            sliceTOTPOP = TOTPOP[:,:,sliceIndex]

            # now can render!
            n = n + 1
            fig1 = plt.figure(1)
            plt.clf()

            def renderActivityFrame(pltObj = plt, actType = activityType):
                # render activity based on activity type
                sliceAct = np.zeros(shape=(xyzSize[d1],xyzSize[d2]), dtype=float)
                if actType == 0: # actSum
                    sliceAct = sliceTOTPOP * particleActivity
                elif actType == 1: # actAv
                    sliceBoolMask = np.isin(sliceTOTPOP, [0], invert=True) # get mask where non-zero is true, zero is false
                    sliceAct = particleActivity * sliceBoolMask
                elif actType == 2: # actSigmoid
                    sliceBoolMask = np.isin(sliceTOTPOP, [0], invert=True) # get mask where non-zero is true, zero is false

                    def sigFnct(val):
                        return 0.5 * (1 - np.tanh((val - actSigAv * (1 + actSigWidth)) / (actSigAv * actSigWidth)))
                    sliceAct = particleActivity * sigFnct(sliceTOTPOP) * sliceBoolMask

                # hack
                if actType == 0:
                    pltObj.imshow(sliceAct,cmap=myMap,vmin=0,aspect=myAspect,origin='lower')
                else: # if using active-av or active-sig, set 1.5 as upper limit
                    pltObj.imshow(sliceAct,cmap=myMap,vmin=0,vmax=1.5,aspect=myAspect,origin='lower')

                #Create the colorbar
                CS3 = None
                if actType == 0:
                    CS3 = pltObj.imshow(sliceAct.T/particleActivity,vmin=0,cmap=myMap,aspect=myAspect)			#pcolor() sucks this is way betterd
                else:
                    CS3 = pltObj.imshow(sliceAct.T/particleActivity,vmin=0,vmax=1.5,cmap=myMap,aspect=myAspect)			#pcolor() sucks this is way betterd
                
                if pltObj != plt:
                    cb=colorbar(CS3, ax=pltObj)
                else:
                    cb=colorbar(CS3)
                if pltObj != plt: # if this is subplots, do a small label
                    cb.ax.set_ylabel(r'$\alpha_\mathrm{C} / \alpha_i$', fontsize = FS)
                else: # otherwise do a full label
                    cb.ax.set_ylabel(r'Activity Strength, $\alpha_\mathrm{C} / \alpha_i$', fontsize = FS)

                # perform matplotlib bits and save
                # title(r'Slice %s=%d'%(sliceDim,sliceIndex), fontsize = FS)
                xlabel(r'$%s$'%labX, fontsize = FS)
                ylabel(r'$%s$'%labY, fontsize = FS)
                plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)

            # if rendering a specific activity type
            if activityType != -1:
                renderActivityFrame()
            else: # if rendering a comparison of different activity levels
                _, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharey=True, sharex=True)

                # density plot
                sliceMax = np.max(sliceTOTPOP) # maximum population in the slice
                ax1.imshow(sliceTOTPOP/sliceMax,cmap=ed.plasma,vmin=0, vmax=1, 
                    aspect=myAspect)
                CS3 = ax1.imshow(sliceTOTPOP.T/sliceMax,vmin=0,vmax=1,cmap=ed.plasma,aspect=myAspect)			#pcolor() sucks this is way better
                cb=colorbar(CS3, ax=ax1)
                cb.ax.set_ylabel(r'$N_C / N_C^\mathrm{max}$', fontsize = FS)
                ax1.set_title(r'Density', fontsize = FS)

                # active sum
                renderActivityFrame(ax2, 0)
                ax2.set_title(r'Sum', fontsize = FS)

                # active av
                renderActivityFrame(ax3, 1)
                ax3.set_title(r'Average', fontsize = FS)

                # active sig
                renderActivityFrame(ax4, 2)
                ax4.set_title(r'Sigmoidal', fontsize = FS)

            name='frame%04d.png'%(n)
            namepdf='frame%04d.pdf'%(n)
            
            ## uncomment below for snapshots!
            plt.axis('off') 
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            if activityType == -1:
                for ax in (ax1, ax2, ax3, ax4):
                    ax.axis('off')
            savefig(name, bbox_inches='tight', dpi=400)
            # save trans pdf
            savefig(namepdf, transparent=True, bbox_inches='tight')

            i=0 # reset tStep line counter
infile.close()

#Animate
print( "\tAnimating ..." )
# single species
name='2DactField_%s%d%s'%(sliceDim,sliceIndex,suffix)
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)

if not keepFrames:
    myCommand="rm s*frame*.png"
    call(myCommand,shell=True)
