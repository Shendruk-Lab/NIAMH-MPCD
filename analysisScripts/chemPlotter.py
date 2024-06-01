import matplotlib
from matplotlib import pyplot as plt
from matplotlib import ticker
from pylab import *
from numpy import ma
from subprocess import call
from numpy import linalg as LA
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors
import os
import csv
import numpy as np
import math
import gyrationTensor_2D
import chemoCurv
###########################################################
### Calculates and plots chemograph of curvature
###########################################################
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
myMap = ed.plasma
FS = 20
test = False
# change following as appropriate - should probs try and make this more sophisticated
repeats = 20
numMono = 24
outT = 200 #will be 4000 later
x = 100 #xdim of system
y = 100
log = False

def myfunc():
    cwd = os.getcwd()
    image3D = np.zeros((repeats, outT, numMono))
    for i in range(1,repeats+1):
        path = str(cwd+"/result"+str(i))
        os.chdir(path)
        chemoCurv.func(path, x, y, numMono,log)
        name = path+'/chemo-curv.dat'
        file = open(name, "r")
        # toss headers
        line = file.readline()
        line = file.readline()
        for t in range(outT):
            for j in range(numMono):
                line = file.readline()
                image3D[i-1][t][j] = float(line.split()[-1])
    # take average along repeats axis
    image2D = np.average(image3D, axis=0)
    print(np.shape(image2D)) #check shape in case averaged along wrong axis
    # std also
    imStd = np.std(image3D, axis=0)
    imErr = imStd/np.sqrt(np.shape(image3D)[0])
    
    # write out
    bend,act, chunks = gyrationTensor_2D.read(cwd)
    name = str(cwd+"/chemCurv_a"+str(act)+"_c"+str(chunks)+".dat")
    outfile = open(name, 'w')
    outfile.write("Average curvature for each time and monomer index with 20 repeats\n")
    outfile.write("time\tindex\tcurv\tstd\tsterr\n")
    for t in range(outT):
        for j in range(numMono):
            outfile.write("%e\t%e\t%e\t%e\t%e\n"%(t*500,j,image2D[t][j],imStd[t][j],imErr[t][j]))
    outfile.close()
    # plot!
    fig,ax = plt.subplots(1)
    plt.cla
    # make sure dir exists
    if log == True:
        figdir = '/home/s1954660/Desktop/summer/analysis/LOGcurv_chemographs'
    else:
        figdir = '/home/s1954660/Desktop/summer/analysis/curv_chemographs'
    figname = str(figdir+"/chemCurv_a"+str(act)+"_c"+str(chunks)+".pdf")
    if log==True:
        chemograph = ax.pcolor(image2D,norm=colors.LogNorm(),cmap=myMap)
    else:
        chemograph = ax.pcolor(image2D,cmap=myMap)
        #chemograph = imshow(image2D,cmap=myMap, origin='lower', aspect='auto')
    cb=fig.colorbar(chemograph)
    if log==True:
        cb.ax.set_ylabel(r"Average Curvature, $\kappa$ - LOGARITHMIC",fontsize=FS)
    else:
        cb.ax.set_ylabel(r"Average Curvature, $\kappa$ (20 repeats)",fontsize=FS)
    xlabel(r"Monomer index",fontsize=FS)
    ylabel(r"time, $t_{\mbox{mpcd}}$/500",fontsize=FS)
    plt.savefig(figname, format='pdf', dpi='figure')
    plt.close(fig)

if test==True:
    # need to start within a chunks folder
    myfunc()
    
else:
    # start in folder cluster01
    clustdir = os.getcwd()
    #datedirlist = os.listdir(clustdir)
    #for date in datedirlist: #go through date folders
    #    datedir = os.path.join(clustdir, date)
        #go through activities
    actdirlist = os.listdir(clustdir)
    for act in actdirlist:
        if 'act' in act:
            actdir = os.path.join(clustdir, act)
            #go through chunks
            chunkdirlist = os.listdir(actdir)
            for chunk in chunkdirlist:
                chunkdir = os.path.join(actdir, chunk)
                os.chdir(chunkdir)
                myfunc()