## ================================================================================
## CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
##
## Copyright (c) 2023 
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be found
## in the LICENSE file at the top level of the distribution and at
## github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
##
## ================================================================================
## Author: Matthew Troemner
## ================================================================================
##
## Description coming soon...
##
##
## ================================================================================

import math
import numpy as np



def overlapCheck(nodes,center,parDiameter,facePoints,binMin,binMax,\
    minPar,maxEdgeLength,aggOffset,parDiameterList):

    # Store particle nodes that fall inside the bin
    binTestParticles = np.all([(nodes[:,0] > binMin[0]) , \
        (nodes[:,0] < binMax[0]) , (nodes[:,1] > binMin[1]) , \
        (nodes[:,1] < binMax[1]) , (nodes[:,2] > binMin[2]) , \
        (nodes[:,2] < binMax[2])],axis=0)
    existingnodes = np.asarray(nodes[binTestParticles,:])
    existingParD = np.asarray(parDiameterList[binTestParticles])

    # Compute distance between particles 
    if len(existingnodes>0):
        nodalDistance = np.linalg.norm(center-existingnodes, axis=1)
        parOffsetDist = nodalDistance - parDiameter/2 - existingParD\
            /2 - aggOffset
    else: 
        parOffsetDist = np.array(([1]))

    # Kill and return if overlap
    if (parOffsetDist<0).any():
        return True,"NA"

    # Store edge nodes that fall inside the bin
    binTestSurf = np.all([(facePoints[:,0] > binMin[0]) , \
        (facePoints[:,0] < binMax[0]) , (facePoints[:,1] > \
            binMin[1]) , (facePoints[:,1] < binMax[1]) ,\
        (facePoints[:,2] > binMin[2]) , (facePoints[:,2] < \
            binMax[2])],axis=0)     
    existingSurf = np.asarray(facePoints[binTestSurf,:])

    # Compute distance between particle and edge nodes
    if len(existingSurf>0):
        surfNodalDistance = np.linalg.norm(center-existingSurf, axis=1)
        parSurfaceDist = np.asarray(surfNodalDistance - parDiameter/2 - 1.1*minPar/2)
    else: 
        parSurfaceDist = np.array(([1]))

    # Kill and return if overlap
    if (parSurfaceDist<0).any():
        return True,"NA"

    # Otherwise return false and check if critically close to surface
    if len(existingSurf>0):
        if (surfNodalDistance<=math.sqrt((maxEdgeLength*math.sqrt(3)/3)**2+(parDiameter/2)**2)).any():
            return False,True

    return False,False