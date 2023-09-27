#plot the avS against time using a _local_ average obtained from directorfield.dat rather than avS.dat
#example usage: 
# python3 ./localS.py 20 ../data1/directorfield.dat data1 ../data2/directorfield.dat data2

import matplotlib.pyplot as plt
import numpy as np
import sys

if(len(sys.argv) < 3):
    print("The script recquires at least 3 arguments. The first should be the warmup time, followed by pairs of the data followed by the label")
    print("Terminating.")
    exit()

#handle arguments - input the path to directorfield.dat for each one you want to render followed by the tag
print( "Arguments:" )
for arg in sys.argv:
	print( "\t" + arg )
STARTUPTIME = float(sys.argv[1]) #how much SIM TIME to wait before starting averaging

#handle graph formatting and style
plt.style.use('shendrukGroupStyle')
import shendrukGroupFormat as ed
MYLW=1.0 #line width

#basic warning if you're passing the wrong format of arguments
if(len(sys.argv) % 2) != 0:
    print("The script expects an ODD number of arguments. The first should be the warmup time, followed by pairs of the data followed by the label")
    print("Terminating.")
    exit()

#### main plotting function ###################################
def plotData(dir: str, lbl: str):
    #load data file
    print("Reading file "+dir+" for data.")
    data = np.loadtxt(dir, delimiter="\t", skiprows=13)
    tData = data[:,0]
    sData = data[:,7]

    #data to plot
    t = []
    avLocalS = []
    #tracking vars
    currTime = 0.0
    sumS = 0.0
    countS = 0
    for i in range(len(tData)):
        if tData[i] < STARTUPTIME: # check for startup time
            continue

        # if the time changed, we need to add this data to our plotting list
        if tData[i] != currTime:
            t.append(currTime)
            avLocalS.append(sumS/countS)

            #reset tracking vars
            currTime = tData[i]
            sumS = 0.0
            countS = 0

        # add data to rolling sum
        countS += 1
        sumS += sData[i]

    plt.plot(t, avLocalS, label = lbl)
###############################################################

N = int((len(sys.argv) / 2) - 1) #number of pairs
print(N)
for i in range(N):
    dataName = sys.argv[i+2]
    dataLabel = sys.argv[i+3]
    plotData(dataName, dataLabel)

plt.xlabel(r"Time, $t$")
plt.ylabel(r"Local Order, $\left\langle S \right\rangle$")
plt.title("Average Local Order Parameter")
plt.legend()
plt.savefig( 'avLovalS.pdf' )
plt.show()