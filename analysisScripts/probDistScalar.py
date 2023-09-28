"""
	Probability density rendering script.

	Created by Tyler Shendruk
"""

from pylab import *
from subprocess import call
from scipy import integrate
import os

###########################################################
### Plots 2D averaging over user defined direction
###########################################################

###########################################################
### Read arguments
###########################################################
print( "Arguments:" )
for arg in sys.argv:
	print( "\t" + arg )
dataName = sys.argv[1]			# Name of the data (should be coarsegrain.dat)
xaxis = sys.argv[2]				# What is this the histogram of?
start = int(sys.argv[3])		# Average after this number
finish = int(sys.argv[4])		# Average before this number
makeMovie = int(sys.argv[5])	# Whether or not to animate temporal data
keepFrames =int(sys.argv[6])	#0=don't keep (delete) frames; 1=keep frames
plotAv = int(sys.argv[7])		# Whether or not to plot average

###########################################################
### Assumed arguments that a user could change
###########################################################
numBins=101
avOutName=dataName.split(".dat")[0]+"_av"

###########################################################
### Format and style
###########################################################
# Use our custom style and colours
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed
#Animation stuff
bitrate=5000
framerate=12		#Number of frames per second in the output video
	# Note that the ideal for this is _very_ variable so play around with it
codec='libx264'		#Other options include mpeg4
suffix='.mp4'

# Figure
fig1 = plt.figure(1)

###########################################################
### Animation
###########################################################
print( 'Reading file for min/max/movie ...' )
if not os.path.isfile(dataName):
	print("%s not found."%dataName)
	exit()
infile = open(dataName,"r")
#toss header
for i in range(13):
	line = infile.readline()
t=0
min=999999.
max=-999999.
maxC=0.0
X=zeros(shape=(numBins),dtype=float)	#Start of the bins
C=zeros(shape=(numBins),dtype=float)
while infile:
	line = infile.readline()
	if( not line ):
		break
	else:
		time=float(line)
		for i in range(numBins):
			line = infile.readline()
			toss,x,c = line.split("\t")
			X[i]=float(x)
			C[i]=float(c)
		line = infile.readline()
		if(t>=start):
			dx=X[1]-X[0]
			min=minimum(min,X[0])
			max=maximum(max,X[-1]+dx)
			maxC=maximum(maxC,C.max())
			if(makeMovie):
				X+=0.5*dx
				fig1 = plt.figure(1)
				plt.clf()
				title(r'Time, %s$\tau$'%(time))
				plot(X,C,'-o')
				xlabel(r'%s'%xaxis)
				ylabel(r'Counts')
				plt.axis(xmax=max, xmin=min, ymax=maxC, ymin=0)
				name='frame%04d.png'%(t)
				savefig(name, bbox_inches='tight')
		t+=1
infile.close()
if(makeMovie):
	#Animate
	print( "\tAnimating ..." )
	name='%s_animation%s'%(avOutName,suffix)
	myCommand="rm %s"%name
	call(myCommand,shell=True)
	myCommand = "ffmpeg -f image2 -r %d"%(framerate)+" -i frame%04d.png"+" -vcodec %s -b %dk -r %d %s"%(codec,bitrate,framerate,name)
	call(myCommand,shell=True)
	if not keepFrames:
		myCommand="rm frame*.png"
		call(myCommand,shell=True)

###########################################################
### Average
###########################################################
print( 'Reading data ...' )
infile = open(dataName,"r")
# Keep header
header=[]
for i in range(13):
	header.append(infile.readline())
t=0
MmM=max-min
avX=linspace(min,max,num=numBins,endpoint=True,dtype=float)	#Start of the bins
avP=zeros(shape=(numBins),dtype=float)
while infile:
	line = infile.readline()
	if( not line ):
		break
	else:
		time=float(line)
		for i in range(numBins):
			line = infile.readline()
			# line = line.split("\t")
			toss,x,c = line.split("\t")
			X[i]=float(x)
			C[i]=float(c)
		line = infile.readline()
		if(t>=start):
			# Integrate under the curve
			dx=X[1]-X[0]
			norm=integrate.simpson(C,dx=dx)
			C/=norm
			X+=0.5*dx
			for i in range(numBins):
				I=int(numBins*(X[i]-min)/MmM)
				avP[I]+=C[i]
		t+=1
infile.close()
dx=avX[1]-avX[0]
norm=integrate.simpson(avP,dx=dx)
avP/=norm

outfile=open(avOutName+".dat", "w")
cutTime=header[-1].split("\t")[1:]
header[-1]="%s\t%s"%(cutTime[0],cutTime[2])
for h in header:
	outfile.write(h)
for i in range(numBins):
	outfile.write("%e\t%e\n"%(avX[i],avP[i]))
outfile.close()

avX+=0.5*dx
if(plotAv):
	fig1 = plt.figure(1)
	plt.clf()
	plot(avX,avP,'-o')
	xlabel(r'%s'%xaxis)
	ylabel(r'PDF')
	plt.axis(xmax=max, xmin=min, ymin=0)
	savefig('%s.png'%(avOutName), bbox_inches='tight')
	savefig('%s.pdf'%(avOutName), bbox_inches='tight')
	show()

exit()