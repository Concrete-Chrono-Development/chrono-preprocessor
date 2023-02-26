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
## Description coming soon...
##
##
## ===========================================================================

import numpy as np



def overlapCheckMPI(center,parDiameter,binMin,binMax,minPar,parOffset,nodes,diameters):

    # Store particle nodes that fall inside a floating bin
    binTestParticles = np.all([(nodes[:,0] > binMin[0]) , \
        (nodes[:,0] < binMax[0]) , (nodes[:,1] > binMin[1]) , \
        (nodes[:,1] < binMax[1]) , (nodes[:,2] > binMin[2]) , \
        (nodes[:,2] < binMax[2])],axis=0)
    existingnodes = nodes[binTestParticles,:]
    existingParD = diameters[binTestParticles]

    # Compute distance between particles 
    if len(existingnodes>0):
        nodalDistance = np.linalg.norm(center-existingnodes, axis=1)
        parOffsetDist = nodalDistance - parDiameter/2 - existingParD\
            /2 - parOffset
    else: 
        parOffsetDist = np.array(([1]))

    # Kill and return if too close to any particle
    if (parOffsetDist<0).any():
        return True

    # Otherwise return false (safe distance from other particles)
    return False