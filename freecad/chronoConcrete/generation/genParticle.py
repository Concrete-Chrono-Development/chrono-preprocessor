


import numpy as np

from freecad.chronoConcrete.generation.particleOverlapCheck     import overlapCheck
from freecad.chronoConcrete.generation.particleInsideCheck     	import insideCheck




def generateParticle(x,facePoints,parDiameter,maxParNum,minC,maxC,\
    vertices,tets,coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,maxPar,\
    aggOffset,verbose,parDiameterList,maxEdgeLength,nodes):
    
    # Generate random numbers to use in generation
    randomN = np.random.rand(newMaxIter*3)
    i=0
    ntet = len(tets)
    # Generate random nodal location
    while True:
        i=i+3

        if i/3 >= newMaxIter:
            i = 0
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

            # Temporarily set this and other instances to check regardless if critical
            # WILL BE FIXED IN FUTURE
            if overlap[1] == True or overlap[1] == False:

                # Check if particle is inside the mesh if critically close          
                inside = insideCheck(vertices,\
                tets,node,parDiameter,binMin,binMax,coord1,coord2,coord3,\
                coord4,maxC)

            else:

                inside = True

            # Indicate placed particle and break While Loop
            if inside == True and overlap[0] == False:
                break


    iterReq = i
    return newMaxIter,node,iterReq