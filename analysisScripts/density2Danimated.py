"""
	Density field rendering script.

	Created by Timofey Kozhukhov
"""

from pylab import *
from subprocess import call

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
normalise=int(sys.argv[11])		#normalise the population (for the cbar) or no
species=int(sys.argv[12])		#1 = do total pop only; 2 = two species
savePDF=int(sys.argv[13])		# 1 for saving transparent pdfs (for papers), 0 for none

###########################################################
### Format and style
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed

# Colour map to use # TODO: adjust this!
myMap=ed.plasma # colour map used when species=1
pop0Col = np.array([0, 1, 0]) # colour to show when ONLY pop0 is present (RGB)
pop1Col = np.array([1, 0, 0]) # colour to show when ONLY pop1 is present (RGB)

def getRGBAColField(totalPop, maxTotalPop, popProp):
	"""
	Helper method to get colours.

	Args:
		totalPop: total number of particles in given cell as a MxN np array
		maxTotalPop: the maximum total number of particles in any cell globally (scalar)
		popProp: proportion of pop0 particles as a MxN np array
	Returns:
		RGBA colour field as a MxNx4 np array
	"""
	# compute alpha
	a = totalPop/maxTotalPop 
	# compute RGB
	r = popProp*pop0Col[0] + (1-popProp)*pop1Col[0]
	g = popProp*pop0Col[1] + (1-popProp)*pop1Col[1]
	b = popProp*pop0Col[2] + (1-popProp)*pop1Col[2]

	return np.stack((r, g, b, a), axis=2)

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

if species < 1:
	print( "Species must be 1 or 2 - not %s"%species )
	exit()
if species > 2:
	print( "More than 2 species are not supported!")
	print( "Species must be 1 or 2 - not %s"%species )
	exit()

# Data
TOTPOP = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
sliceTOTPOP = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)

POP0, slicePop0, POP1, slicePOP1, slicePOP0prop, sliceRGBA = None, None, None, None, None, None
if species == 2: # lets be memory efficient
	POP0 = zeros(shape=(xyzSize[0],xyzSize[1],xyzSize[2]),dtype=float)
	slicePOP0 = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float)
	slicePOP0prop = zeros(shape=(xyzSize[d1],xyzSize[d2]),dtype=float) # proportion of population 0 (POP0/TOTPOP)
	sliceRGBA = zeros(shape=(xyzSize[d1],xyzSize[d2],4),dtype=float) # RGBA colour for each cell

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

			if species == 2: # now handle multispecies info
				sp0 = currLine[8]
				POP0[int(Qx)][int(Qy)][int(Qz)] = float(sp0) # this is all we need for now

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
			sliceMax = np.max(sliceTOTPOP) # maximum population in the slice
			if species == 2: # for multispecies
				slicePop0 = POP0[:,:,sliceIndex]
				slicePop0prop = slicePop0/sliceTOTPOP

			# now can render!
			n = n + 1
			fig1 = plt.figure(1)
			plt.clf()

			# handle normalisation if required
			renderTOTPOP = sliceTOTPOP
			if normalise:
				renderTOTPOP = renderTOTPOP/sliceMax

			#handle single species render first
			imshow(renderTOTPOP,cmap=myMap,vmin=0, aspect=myAspect)

			#Create the colorbar
			CS3 = imshow(renderTOTPOP.T,vmin=0,cmap=myMap,aspect=myAspect)			#pcolor() sucks this is way better
			cb=colorbar(CS3)
			if normalise:
				cb.ax.set_ylabel(r'$N_C / N_C^\mathrm{max}$', fontsize = FS)
			else:
				cb.ax.set_ylabel(r'$N_C$', fontsize = FS)

			# perform matplotlib bits and save
			# title(r'Slice %s=%d'%(sliceDim,sliceIndex), fontsize = FS)
			xlabel(r'$%s$'%labX, fontsize = FS)
			ylabel(r'$%s$'%labY, fontsize = FS)
			plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)

			name='s=1frame%04d.png'%(n)
			namepdf='s=1frame%04d.pdf'%(n)
			
			## uncomment below for snapshots!
			plt.axis('off') 
			plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
			savefig(name, bbox_inches='tight')
			# save trans pdf
			if savePDF: savefig(namepdf, transparent=True, bbox_inches='tight')

			# savefig( name )
				
			# also handle multispecies render if necessary
			if species == 2: # for multispecies
				plt.clf()
				sliceRGBA = getRGBAColField(sliceTOTPOP, sliceMax, slicePop0prop)
				imshow(sliceRGBA, aspect=myAspect)
				# perform matplotlib bits and save
				# title(r'Slice %s=%d'%(sliceDim,sliceIndex), fontsize = FS)
				xlabel(r'$%s$'%labX, fontsize = FS)
				ylabel(r'$%s$'%labY, fontsize = FS)
				plt.axis(xmax=xyzSize[d1], xmin=0, ymax=xyzSize[d2], ymin=0)
				name='s=2frame%04d.png'%(n)
				savefig( name )
			

			i=0 # reset tStep line counter
infile.close()

#Animate
print( "\tAnimating ..." )
# single species
name='2Ddens_%s%d%s'%(sliceDim,sliceIndex,suffix)
myCommand="rm %s"%name
call(myCommand,shell=True)
myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i s=1frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
call(myCommand,shell=True)

if species == 2: # for multispecies
	name='2Ddens_s=%d_%s%d%s'%(species,sliceDim,sliceIndex,suffix)
	myCommand="rm %s"%name
	call(myCommand,shell=True)
	myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i s=2frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
	call(myCommand,shell=True)

if not keepFrames:
	myCommand="rm s*frame*.png"
	call(myCommand,shell=True)
