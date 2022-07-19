# A class to handle defects for the purposes of drawing

from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt

@dataclass
class Defect:
    pos : np.ndarray
    charge : float
    angle : float
    defectGroup : int = -1

    def drawDefect(self, lineLength : float, lineWidth : float, stylePlus="r", 
                    styleMinus="b", pltObj = plt):
        """
        Draws a defect on a plot.

        :param lineLength: length of the lines used to draw the defect
        :type lineLength: float
        :param lineWidth: width of the lines used to draw the defect
        :type lineWidth: float
        :param stylePlus: style of lines used to draw +1/2 defects
        :param styleMinus: style of lines used to draw -1/2 defects
        :param pltObj: The axes to draw on
        """
        ##FIXME: only works for 2d
        def drawLine(pos, angleRelDefect, length, style):
            x = pos[0] + length*np.cos(angleRelDefect+self.angle)
            y = pos[1] + length*np.sin(angleRelDefect+self.angle)
            return pltObj.plot([pos[0], x], [pos[1], y], style, lw=lineWidth)

        if self.charge > 0:
            # +1/2 defect
            drawLine(self.pos-0.5*lineLength*np.array(
                [np.cos(self.angle), np.sin(self.angle), 0]), 
                np.pi/6, lineLength * 1.5, stylePlus)
            drawLine(self.pos-0.5*lineLength*np.array(
                [np.cos(self.angle), np.sin(self.angle), 0]), 
                -np.pi/6, lineLength * 1.5, stylePlus)
        elif self.charge < 0: 
            # -1/2 defect
            drawLine(self.pos, 0, lineLength, styleMinus)
            drawLine(self.pos, 2*np.pi/3, lineLength, styleMinus)
            drawLine(self.pos, 4*np.pi/3, lineLength, styleMinus)

    def distTo(self, otherDefect, l):
        """
        Get distance to another defect

        :param otherDefect: Defect to get distance to
        :param l: box dimensions, can be 3d array
        """
        diff = self.pos - otherDefect.pos
        return np.mod(diff + 0.5* l, l) - 0.5*l # MIC

def extractData(dir: str, sysDim):
    """
    Extracts defect data from a topology data file.

    :param dir: Path to the topology data file
    :type dir: str
    :param sysDim: Dimensions of the system domain
    """

    #load data file, cannot use numpy for this because of the data format
    file = dir
    infile = open(file, "r")
    for i in range(13): #toss header
        line = infile.readline()

    #data we're going to plot
    TDATA = [0]
    DEFECTLIST = [] # list of defects on every timestep

    lastT = 0. #last Time measured
    DEFECTDATA = [] #list of defects on the current timestep
    
    DEFECTIDARRAY = -1*np.ones((sysDim[0], sysDim[1], sysDim[2])) 
    # 3D array containing defect ID info, [x][y][z], default -1

    error = False # error catching
    currLine = 14
    while infile: #loop through file to look for the target data we want
        line = infile.readline()

        def endTStep(): # code to be run at the end of each timestep
            collatedDefect = collapseDefects(DEFECTDATA, DEFECTIDARRAY, 
                sysDim)
            DEFECTLIST.append(collatedDefect)

        if (not line): #leave loop if EoF
            endTStep()
            break
        else:
            try:
                t,qx,qy,qz,charge,angle = line.split("\t", 6)
                t = float(t)
                charge = float(charge)
            except:
                error = True
                print("Error reading file "+dir+" on line +"+str(currLine))
                continue

            if t > lastT: # if there's a new timestep worth of data, prepare for the next block of data
                if not error: endTStep() # only append data if not error
                if not error: TDATA.append(t) 
                lastT = t                
                error = False

                DEFECTDATA = [] #clear defect data list
                DEFECTIDARRAY = -1*np.ones((sysDim[0], sysDim[1], sysDim[2]))  
                # 3D array containing defect ID info, [x][y][z], default -1

            if abs(abs(charge)-0.5) < 0.001: # if we have charge \pm 1/2 then we have a defect
                DEFECTDATA.append(Defect(pos=np.array([float(qx), float(qy), 
                    float(qz)]), charge=float(charge), angle=float(angle)))

                #add to the defect ID array too
                DEFECTIDARRAY[int(qx)][int(qy)][int(qz)] = len(DEFECTDATA) - 1
        
        currLine += 1
    
    return TDATA, DEFECTLIST

#merge neighbouring defects of like charge together
def collapseDefects(DEFECTDATA, DEFECTIDARRAY, sysDim):
    # a group of defects to be collated
    class DefectGroup:
        def __init__(self, groupID):
            self.charge = 0
            self.posList = [[], [], []] # stores as [index][elem]
            self.angles = []
            self.ID = groupID

        def addDefect(self, target : Defect):
            if self.charge != 0: # if charge is already set, do a check
                if self.charge != target.charge:
                    print("Trying to assign a defect of incorrect charge to group ID "+self.ID)
            else: #if charge not set, set charge
                self.charge = target.charge

            # add to posList and increment count
            for i,coord in enumerate(target.pos):
                self.posList[i].append(coord)

            # add to angle list
            self.angles.append(target.angle)

            # tag defect
            target.defectGroup = self.ID

        def getDefect(self): #convert a defect group back to a defect
            # get av x, y, z pos
            avX = sum(self.posList[0])/len(self.posList[0])
            avY = sum(self.posList[1])/len(self.posList[1])
            avZ = sum(self.posList[2])/len(self.posList[2])

            # get av angle - explicit average wont work in all cases
            unitVecList = []
            for angle in self.angles:
                unitVecList.append(np.array([np.cos(angle), np.sin(angle)]))
            avUnitVec = sum(unitVecList)/len(unitVecList)
            avAng = np.arctan2(avUnitVec[1], avUnitVec[0])

            return Defect(pos=np.array([avX, avY, avZ]), charge=self.charge, 
                angle=avAng)

    currGroupCount = 0 #init number of groups
    groupList = [] # list of current groups

    #loop through defects
    for i, defect in enumerate(DEFECTDATA):
        # method to check (and potentially add) neighbour to groupList
        def checkNeighbour(coord):
            # check list for neighbour ID
            qx = int(coord[0])
            qy = int(coord[1])
            qz = int(coord[2])
            nID = int(DEFECTIDARRAY[qx][qy][qz] )

            if nID != -1: #if there exists a defect there...
                if DEFECTDATA[nID].charge == defect.charge: #... of like charge
                    nDGroup = DEFECTDATA[nID].defectGroup #get defect group
                    
                    if nDGroup != -1: #if nDGroup is assigned
                        groupList[nDGroup].addDefect(defect) #assign defect

                        return True # return true
                    else: return False
                else: return False
            else: return False #if no defect is there, return false

        # construct coordDiffs to check
        coord = defect.pos
        coordDiffs = []
        for xd in range(-1, 2):
            for yd in range(-1, 2):
                for zd in range(-1, 2):
                    x = coord[0]+xd
                    y = coord[1]+yd
                    z = coord[2]+zd

                    # check out of bounds
                    if (x < 0) or (x >= sysDim[0]):
                        continue
                    if (y < 0) or (y >= sysDim[1]):
                        continue
                    if (z < 0) or (z >= sysDim[2]):
                        continue

                    coordDiffs.append([x, y, z]) # if here then add
                
        status = False # whether we were assigned or not
        for coord in coordDiffs:
            status = checkNeighbour(coord) # consider each neighbours
            if status == True: # if succesfull then move onto next defect
                break
        
        #if we're here and still not assigned, then create a new group
        if status == False:
            groupList.append(DefectGroup(currGroupCount))
            groupList[-1].addDefect(defect)

            currGroupCount += 1 #increment group counter

    #export a dumped list of defects
    return [dg.getDefect() for dg in groupList]