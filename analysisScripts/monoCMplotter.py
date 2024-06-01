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
import monoContactMap
# Use our custom style
plt.style.use('shendrukGroupStyle')
# Use our custom colours
import shendrukGroupFormat as ed
# Colour map to use
myMap = ed.viridis
FS = 20
test = False
log = True

if test==True:
    # need to be within a chunks folder
    bend, act, chunks = gyrationTensor_2D.read(os.getcwd())
    monoContactMap.func(os.getcwd(), bend, act, chunks, log)
else:
    # start within cluster01 (for example)
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
                bend, act, chunks = gyrationTensor_2D.read(chunkdir)
                monoContactMap.func(chunkdir, bend, act, chunks, log)