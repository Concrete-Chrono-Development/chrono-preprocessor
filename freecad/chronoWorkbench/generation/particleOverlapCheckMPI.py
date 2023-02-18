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

import numpy as np



def overlapCheckMPI(center,aggDiameter,binMin,binMax,minPar,aggOffset,nodes,diameters):

    # Store particle nodes that fall inside the bin
    binTestParticles = np.all([(nodes[:,0] > binMin[0]) , \
        (nodes[:,0] < binMax[0]) , (nodes[:,1] > binMin[1]) , \
        (nodes[:,1] < binMax[1]) , (nodes[:,2] > binMin[2]) , \
        (nodes[:,2] < binMax[2])],axis=0)
    existingnodes = nodes[binTestParticles,:]
    existingAggD = diameters[binTestParticles]

    # Compute distance between particles 
    if len(existingnodes>0):
        nodalDistance = np.linalg.norm(center-existingnodes, axis=1)
        aggOffsetDist = nodalDistance - aggDiameter/2 - existingAggD\
            /2 - aggOffset
    else: 
        aggOffsetDist = np.array(([1]))

    # Kill and return if overlap
    if (aggOffsetDist<0).any():
        return True

    # Otherwise return false
    return False