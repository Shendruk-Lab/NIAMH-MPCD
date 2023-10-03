import matplotlib.colors as mcolors

def colourListRGB(  ):
	return [ [146./256.,26./256.,28./256.],[205./256.,55./256.,39./256.],[225./256.,175./256.,37./256.],[108./256.,181./256.,156./256.],[107./256.,107./256.,177./256.],[109./256.,110./256.,113./256.],[35./256.,31./256.,32./256.] ]

def invColourListRGB(  ):
	c=colourListRGB(  )
	return c[::-1]

def colourList4Map(clist):
	colourList=[]
	for i in range(len(clist)-1):
		colourList.append(clist[i])
		colourList.append(clist[i+1])
		if i!=len(clist)-2:
			colourList.append(float(i+1)/float(len(clist)-1))
	return colourList

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def make_colormap_withWhite(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    cdict['red'][0][2]=1.0
    cdict['green'][0][2]=1.0
    cdict['blue'][0][2]=1.0
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)
