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


def overlapCheck(nodes, center, parDiameter, facePoints, binMin, binMax,
    minPar, maxEdgeLength, aggOffset, parDiameterList):

    # Combine the conditions for all particles and store the result in an array
    binTestParticles = np.logical_and.reduce([(nodes[:,0] > binMin[0]),
        (nodes[:,0] < binMax[0]), (nodes[:,1] > binMin[1]),
        (nodes[:,1] < binMax[1]), (nodes[:,2] > binMin[2]),
        (nodes[:,2] < binMax[2])])

    # Store particle nodes that fall inside the bin
    existingNodes = nodes[binTestParticles, :]
    existingParD = parDiameterList[binTestParticles]

    # Combine the conditions for all surface points and store the result in an array
    binTestSurf = np.all([(facePoints[:,0] > binMin[0]),
        (facePoints[:,0] < binMax[0]), (facePoints[:,1] > binMin[1]),
        (facePoints[:,1] < binMax[1]), (facePoints[:,2] > binMin[2]),
        (facePoints[:,2] < binMax[2])], axis=0)

    # Store edge nodes that fall inside the bin
    existingSurf = facePoints[binTestSurf, :]

    if existingNodes.shape[0] > 0:
        nodalDistance = np.linalg.norm(center - existingNodes, axis=1)
        parOffsetDist = nodalDistance - parDiameter/2 - existingParD/2 - aggOffset
        if (parOffsetDist < 0).any():
            return True, "NA"
    else:
        parOffsetDist = np.array([1])

    if existingSurf.shape[0] > 0:
        surfNodalDistance = np.linalg.norm(center - existingSurf, axis=1)
        if (surfNodalDistance**2 <= (maxEdgeLength * np.sqrt(3) / 3)**2 + (parDiameter / 2)**2).any():
            return False, True
        if (surfNodalDistance - parDiameter/2 - 1.1 * minPar / 2 < 0).any():
            return True, "NA"
    else:
        parSurfaceDist = np.array([1])

    return False, False
