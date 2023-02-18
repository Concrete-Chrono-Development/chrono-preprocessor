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

from freecad.chronoWorkbench.generation.particleOverlapCheck     import overlapCheck
from freecad.chronoWorkbench.generation.particleInsideCheck     	import insideCheck


def generateParticleMPI(facePoints,maxParNum,minC,maxC,\
    vertices,tets,coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,maxPar,\
    aggOffset,verbose,parDiameterList,maxEdgeLength,max_dist,nodes,parDiameter):
    
    # Generate random numbers to use in generation
    randomN = np.random.rand(newMaxIter*3)
    i=0
    ntet = len(tets)
    # Generate random nodal location
    while True:
        i=i+3

        if i/3 >= newMaxIter:
            i = 3
            newMaxIter = newMaxIter*2
            randomN = np.random.rand(newMaxIter*3)

        if newMaxIter >= maxIter:
            print("This particle has exceeeded the %r specified maximum iterations allowed." % (maxIter))
            print('Now exitting...')
            exit()

        # Random point selection in random tet prism container
        tetVerts = np.vstack((vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),0]-1),:],\
            vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),1]-1),:],\
            vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),2]-1),:],\
            vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),3]-1),:]))

        tetMin = np.amin(tetVerts, axis=0)
        tetMax = np.amax(tetVerts, axis=0)

        node = np.array([randomN[i]*(tetMax[0]-tetMin[0])+tetMin[0],\
            randomN[i+1]*(tetMax[1]-tetMin[1])+tetMin[1],randomN[i+2]\
            *(tetMax[2]-tetMin[2])+tetMin[2]]).T
        node = node[np.newaxis,:]           


        # Obtain extents for floating bin
        binMin = np.array(([node[0,0]-parDiameter/2-maxPar/2-aggOffset,\
            node[0,1]-parDiameter/2-maxPar/2-aggOffset,node[0,2]-\
            parDiameter/2-maxPar/2-aggOffset]))
        binMax = np.array(([node[0,0]+parDiameter/2+maxPar/2+aggOffset,\
            node[0,1]+parDiameter/2+maxPar/2+aggOffset,node[0,2]+\
            parDiameter/2+maxPar/2+aggOffset]))


        # Check if particle overlapping any existing particles or bad nodes
        overlap = overlapCheck(nodes,node,parDiameter,facePoints,binMin,\
            binMax,minPar,maxEdgeLength,aggOffset,parDiameterList)

        if overlap[0] == False:
            
            # If critically close to the surface:
            if overlap[1] == True:


                # Check if particle is inside the mesh if critically close          
                inside = insideCheck(vertices,tets,node,parDiameter,binMin,binMax,coord1,\
                                    coord2,coord3,coord4)

            else:

                inside = True

            # Indicate placed particle and break While Loop
            if inside == True and overlap[0] == False:
                break


    return np.append(node[0,:],[parDiameter,newMaxIter,int(i/3)])