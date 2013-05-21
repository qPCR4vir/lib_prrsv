#!/usr/local/bin/python 
#
# mpl.py
#
# matlplotlib snippets
#

__verison__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

def c1Kc2CM(color1,color2,midPoint=0.5,segments=256):
    """returns a color map like

    0.0 --- midpoint ---- 1.0
    color1-> black  -> color2

    where colors are RGB tuples with values 0<->1.
    """

    s,m,e= (0.0, float(midpoint), 1.0)
    #TODO: finish




def greenRedCM():
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (0.5, 0.0, 0.0),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 1.0, 1.0),
                       (0.5, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.0, 0.0),
                      (0.5, 0.0, 0.0),
                      (1.0, 0.0, 0.0))}
    return matplotlib.colors.LinearSegmentedColormap('green2red',cdict,256)

def redGrayBlueCM():
    cdict = {'red': ((0.0, 0.002, 0.002),
                     (0.45, 0.3, 0.3),
                     (0.5, 0.2, 0.2),
                     (0.55, 0.3, 0.3), 
                     (0.75, 0.85, 0.85),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.009, 0.009),
                       (0.25, 0.0, 0.0),
                       (0.45, 0.3, 0.3),
                       (0.5, 0.2, 0.2),
                       (0.55, 0.3, 0.3),              
                       (0.75, 0.65, 0.65),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.65, 0.65),
                      (0.25, 0.8, 0.8),
                      (0.45, 0.3, 0.3),
                      (0.5, 0.2, 0.2),
                      (0.55, 0.3, 0.3),
                      (0.75, 0.02, 0.02),
                      (1.0, 0.0, 0.0))}
    return matplotlib.colors.LinearSegmentedColormap('red2gray2blue',
                                                     cdict, 256)
