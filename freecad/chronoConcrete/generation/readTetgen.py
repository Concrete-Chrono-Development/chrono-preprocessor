import numpy as np


def readTetgen(nodeFile, tetFile):                                       

    allNodes = np.loadtxt(nodeFile, usecols=(1,2,3), \
        skiprows=1)                                   

    allTets = np.loadtxt(tetFile, usecols=(1,2,3,4), \
        skiprows=1)   

    return allNodes, allTets