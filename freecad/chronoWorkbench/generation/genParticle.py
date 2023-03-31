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
## This file contains the function to generate a particle and outputs the
## location of the particle as well as the maximum number of iterations
## allowed and the number of iterations required to place the particle.
##
## ===========================================================================

import numpy as np

from freecad.chronoWorkbench.generation.checkParticleOverlap    import checkParticleOverlap
from freecad.chronoWorkbench.generation.checkParticleInside     import checkParticleInside


def generateParticle(facePoints,parDiameter,\
    vertices,tets,newMaxIter,maxIter,minPar,maxPar,\
    parOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,nodes):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - facePoints:       List of points on the surface of the mesh
    - parDiameter:      Diameter of the particle
    - vertices:         List of vertices of the mesh
    - tets:             List of tetrahedrons of the mesh
    - newMaxIter:       Maximum number of iterations to try to place a particle
    - maxIter:          Maximum number of iterations to try to place a particle
    - minPar:           Minimum particle diameter
    - maxPar:           Maximum particle diameter
    - parOffset:        Offset coefficient for particle placement
    - parDiameterList:  List of particle diameters
    - coord1:           Coordinate 1 of the tets
    - coord2:           Coordinate 2 of the tets
    - coord3:           Coordinate 3 of the tets
    - coord4:           Coordinate 4 of the tets
    - maxEdgeLength:    Maximum edge length of the mesh
    - max_dist:         Maximum distance from the surface
    - nodes:            List of nodes
    --------------------------------------------------------------------------
    ### Outputs ###
    - node:             Node location of the particle
    - newMaxIter:       Maximum number of iterations to try to place a particle
    - iterReq:          Number of iterations required to place a particle
    --------------------------------------------------------------------------
    """  

    # Generate random numbers to use in generation
    randomN = np.random.rand(newMaxIter*3)    

    # Generate random nodal location
    iterReq = 0
    while True:
        iterReq = iterReq + 3

        if iterReq/3 >= newMaxIter:
            iterReq = 0
            newMaxIter = newMaxIter * 2
            randomN = np.random.rand(newMaxIter*3)

        if newMaxIter >= maxIter:
            print("This particle has exceeeded the %r specified maximum iterations allowed." % (maxIter))
            print('Now exitting...')
            exit()

        # Random point selection in random tet prism container    
        tetIndex = int(np.around(randomN[iterReq] * len(tets))) - 1
        tetVerts = vertices[tets[tetIndex]-1]

        tetMin = np.amin(tetVerts, axis=0)
        tetMax = np.amax(tetVerts, axis=0)

        node = randomN[iterReq:iterReq+3] * (tetMax - tetMin) + tetMin
        node = node[np.newaxis,:]

        # Obtain extents for floating bin
        binMin = node[0,:] - parDiameter/2 - maxPar/2 - parOffset
        binMax = node[0,:] + parDiameter/2 + maxPar/2 + parOffset

        # Check if particle overlapping any existing particles or bad nodes
        overlap = checkParticleOverlap(nodes,node,parDiameter,facePoints,binMin,\
            binMax,minPar,maxEdgeLength,parOffset,parDiameterList)

        # If does not overlap an existing particle set overlap[0] = False
        if overlap[0] == False:
            
            # If critically close to the surface set overlap[1] = True
            if overlap[1] == True:

                # Check if particle is inside the mesh if critically close          
                inside = checkParticleInside(vertices,tets,node,parDiameter,binMin,binMax,coord1,\
                                    coord2,coord3,coord4)

            else:
                inside = True

            # Indicate placed particle and break While Loop
            if inside == True and overlap[0] == False:
                return newMaxIter,node,iterReq
