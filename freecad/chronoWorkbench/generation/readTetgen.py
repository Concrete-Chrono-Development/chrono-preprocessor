## ===========================================================================
## CHRONO WORKBENCH:github.com/Concrete-Chrono-Development/chrono-preprocessor
##
## Copyright (c) 2023 
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be
## found in the LICENSE file at the top level of the distribution and at
## github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
##
## ===========================================================================
## Developed by Northwestern University
## Primary Authors: Matthew Troemner
## ===========================================================================
##
## This functions reads the output node and tet file from Tetgen and 
## constructs NumPy arrays with all of the nodes and tetrahedra.
##
## ===========================================================================

import numpy as np


def readTetgen(nodeFile, tetFile, edgeFile):                                       

    """
    Variable List:
    --------------------------------------------------------------------------
    ### Inputs ###
    nodeFile:        file path of the node file to read
    tetFile:         file path of the tetrahedron file to read
    --------------------------------------------------------------------------
    ### Outputs ###
    allNodes:        numpy array with node coordinates for all tetrahedra
    allTets:         numpy array with the vertex indices of each tetrahedron
    --------------------------------------------------------------------------
    """

    allNodes = np.loadtxt(nodeFile, usecols=(1,2,3),skiprows=1)                                   
    allTets = np.loadtxt(tetFile, usecols=(1,2,3,4),skiprows=1)  
    allEdges = np.loadtxt(edgeFile, usecols=(1,2),skiprows=1) 

    return allNodes, allTets, allEdges