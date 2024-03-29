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
## For U.S. Army ERDC Contract No. W9132T22C0015
## Primary Authors: Matthew Troemner
## ===========================================================================
##
## This file contains the function to check if a point is inside a tet
##
## ===========================================================================

import numpy as np

def check_LDPMCSL_pointInside(vertices,tets,center,binMin,binMax,coord1,\
    coord2,coord3,coord4):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - vertices:         Nodes of the tets
    - tets:             Tets of the mesh
    - center:           Point coordinates
    - binMin:           Minimum coordinates of the bin
    - binMax:           Maximum coordinates of the bin
    - coord1:           Coordinates of the first vertex of each tet
    - coord2:           Coordinates of the second vertex of each tet
    - coord3:           Coordinates of the third vertex of each tet
    - coord4:           Coordinates of the fourth vertex of each tet
    --------------------------------------------------------------------------
    ### Outputs ###
    - Boolean:          True if the point is inside tets, False if not
    --------------------------------------------------------------------------
    """  

    # Convert center to a 1D array
    center = center.flatten()

    # Store tet vertices that fall inside the bin
    coord1 = np.all([(coord1[:,0] > binMin[0]) , (coord1[:,0] < binMax[0]),\
                     (coord1[:,1] > binMin[1]) , (coord1[:,1] < binMax[1]),\
                     (coord1[:,2] > binMin[2]) , (coord1[:,2] < binMax[2])],axis=0)      
    coord2 = np.all([(coord2[:,0] > binMin[0]) , (coord2[:,0] < binMax[0]),\
                     (coord2[:,1] > binMin[1]) , (coord2[:,1] < binMax[1]),\
                     (coord2[:,2] > binMin[2]) , (coord2[:,2] < binMax[2])],axis=0)          
    coord3 = np.all([(coord3[:,0] > binMin[0]) , (coord3[:,0] < binMax[0]),\
                     (coord3[:,1] > binMin[1]) , (coord3[:,1] < binMax[1]),\
                     (coord3[:,2] > binMin[2]) , (coord3[:,2] < binMax[2])],axis=0)  
    coord4 = np.all([(coord4[:,0] > binMin[0]) , (coord4[:,0] < binMax[0]),\
                     (coord4[:,1] > binMin[1]) , (coord4[:,1] < binMax[1]),\
                     (coord4[:,2] > binMin[2]) , (coord4[:,2] < binMax[2])],axis=0)  

    # Get tets with a vertex that falls inside the bin
    binTets = np.any([coord1,coord2,coord3,coord4],axis=0)

    # Store vertices of those tets
    coord1 = vertices[tets.astype(int)[binTets,0]-1]
    coord2 = vertices[tets.astype(int)[binTets,1]-1]
    coord3 = vertices[tets.astype(int)[binTets,2]-1]
    coord4 = vertices[tets.astype(int)[binTets,3]-1]


    # Check whether the point falls inside of tets
    inside = 0
    emptyOnes = np.ones(len(coord1[:,0]))

    D00 = np.rot90(np.dstack((coord1[:,0],coord1[:,1],coord1[:,2],\
        emptyOnes)), 3)
    D01 = np.rot90(np.dstack((coord2[:,0],coord2[:,1],coord2[:,2],\
        emptyOnes)), 3)
    D02 = np.rot90(np.dstack((coord3[:,0],coord3[:,1],coord3[:,2],\
        emptyOnes)), 3)
    D03 = np.rot90(np.dstack((coord4[:,0],coord4[:,1],coord4[:,2],\
        emptyOnes)), 3)

    D0 = np.linalg.det(np.hstack((D00,D01,D02,D03)))

    D10 = np.rot90(np.dstack((emptyOnes*center[0],\
        emptyOnes*center[1],emptyOnes*center[2],emptyOnes)), 3)
    
    D1 = np.linalg.det(np.hstack((D10,D01,D02,D03)))
    D2 = np.linalg.det(np.hstack((D00,D10,D02,D03)))
    D3 = np.linalg.det(np.hstack((D00,D01,D10,D03)))
    D4 = np.linalg.det(np.hstack((D00,D01,D02,D10)))

    if np.logical_and(np.logical_and(np.sign(D0) == np.sign(D1),\
        np.sign(D0) == np.sign(D2)),\
        np.logical_and(np.sign(D0) == np.sign(D3),\
        np.sign(D0) == np.sign(D4))).any():
        return True
    
    else:
        return False
