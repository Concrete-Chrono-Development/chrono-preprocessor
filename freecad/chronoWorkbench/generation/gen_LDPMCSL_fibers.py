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
## 
## ===========================================================================
## 
## This file contains the function to generate a fiber and outputs the
## location of the fiber as well other fiber properties. 
##
## ===========================================================================

import numpy as np



def gen_LDPMCSL_fibers(vertices,tets,coord1,coord2,coord3,coord4,maxIter,\
    lFiber,maxC,maxAggD,fiberOrientation,orientationStrength,triangles,\
    cutFiber):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - vertices:             List of vertices of the mesh
    - tets:                 List of tetrahedrons of the mesh
    - coord1:               Coordinate 1 of the tets
    - coord2:               Coordinate 2 of the tets
    - coord3:               Coordinate 3 of the tets
    - coord4:               Coordinate 4 of the tets
    - maxIter:              Maximum number of iterations to try to place a fiber
    - lFiber:               Length of the fiber
    - maxC:                 Maximum distance from the surface
    - maxAggD:              Maximum aggregate diameter
    - fiberOrientation:     Orientation of the fiber
    - orientationStrength:  Strength of the orientation
    - triangles:            List of triangles of the mesh
    - cutFiber:             Option to cut the fiber
    --------------------------------------------------------------------------
    ### Outputs ###
    - p1Fiber:              First node of the fiber
    - p2Fiber:              Second node of the fiber
    - orienFiber:           Orientation of the fiber
    - lFiber:               Length of the fiber
    --------------------------------------------------------------------------
    """



    # Generate random numbers to use in generation
    randomN1 = np.random.rand(maxIter*3)
    randomN2 = np.random.rand(maxIter*3)
    iterReq=0

    # Generate random nodal location
    while True:

        iterReq = iterReq + 3

        if iterReq >= len(randomN1):
            print("This fiber has exceeeded the %r specified maximum iterations allowed." % (maxIter))
            print('Now exitting...')
            exit()

        # Random point selection in random tet prism container
        tetIndex = int(np.around(randomN1[iterReq] * len(tets))) - 1
        tetVerts = vertices[tets[tetIndex]-1]

        tetMin = np.amin(tetVerts, axis=0)
        tetMax = np.amax(tetVerts, axis=0)

        p1Fiber = np.array([randomN1[iterReq]*(tetMax[0]-tetMin[0])+tetMin[0],\
            randomN1[iterReq+1]*(tetMax[1]-tetMin[1])+tetMin[1],randomN1[iterReq+2]\
            *(tetMax[2]-tetMin[2])+tetMin[2]]).T        

        if fiberOrientation == []:

            # Option for Totally Random Orientation (Get spherical -> Cartesian -> Normalize)

            orienFiber1 = np.array((1,randomN2[iterReq+1]*2*np.pi,randomN2[iterReq+2]*np.pi))

            orienFiber2 = np.array((np.sin(orienFiber1[2])*np.cos(orienFiber1[1]),np.sin(orienFiber1[2])*np.sin(orienFiber1[1]),np.cos(orienFiber1[2])))

            orienFiber = np.array((orienFiber2[0],orienFiber2[1],orienFiber2[2]))/\
                np.linalg.norm(np.array((orienFiber2[0],orienFiber2[1]\
                ,orienFiber2[2])))

        else:

            # Option with Preferred Orientation

            strength = (6**(4-4*orientationStrength)-1)/200
        
            v = np.empty(2)
            j = 0

            while j < 2:
                y = np.random.normal(0, strength, 1)
                if y > -1 and y < 1:
                    v[j] = y
                    j = j+1

            # Normalize fiber orientation
            orienFiber1 = np.array((fiberOrientation[0],fiberOrientation[1],fiberOrientation[2]))/\
                np.linalg.norm(np.array((fiberOrientation[0],fiberOrientation[1],fiberOrientation[2])))

            # Get spherical coordinates
            orienFiber2 = np.array((1,np.arctan2(orienFiber1[1],orienFiber1[0]),np.arccos(orienFiber1[2]/(orienFiber1[0]**2+orienFiber1[1]**2+orienFiber1[2]**2)**0.5)))

            # Perturb values
            orienFiber3 = np.array((1,orienFiber2[1]+np.pi*v[0],orienFiber2[2]+np.pi/2*v[1]))

            # Convert back to Cartesian
            orienFiber = np.array((np.sin(orienFiber3[2])*np.cos(orienFiber3[1]),np.sin(orienFiber3[2])*np.sin(orienFiber3[1]),np.cos(orienFiber3[2])))

            randSign = np.random.rand(1)
            if randSign<0.5:
                sign = -1
            else:
                sign = 1

            # Include opposite direction
            orienFiber = sign*orienFiber


        p2Fiber = p1Fiber+orienFiber*lFiber

        # Obtain extents for floating bin
        binMin = np.amin(np.vstack((p1Fiber,p2Fiber)), axis=0)-maxAggD-lFiber
        binMax = np.amax(np.vstack((p1Fiber,p2Fiber)), axis=0)+maxAggD+lFiber

        # Check if fiber is inside the mesh    
        inside = False     
        inside = insideCheckFiber(vertices,tets,p1Fiber,p2Fiber,\
            binMin,binMax,coord1,coord2,coord3,coord4,maxC)

        # Indicate placed fiber and break While Loop
        if inside == True:
            return p1Fiber, p2Fiber, orienFiber, lFiber
        
        # Find point fiber intersects external surface and trim accordingly
        else:

            # Find point fiber intersects external surface and trim accordingly
            if cutFiber in ['on','On','Y','y','Yes','yes']:                  

                # Get all surface triangle coordinates
                triangles = triangles.astype(int)
            
                coords0 = vertices[triangles[:,0]-1]
                coords1 = vertices[triangles[:,1]-1]
                coords2 = vertices[triangles[:,2]-1] 

                averageTriangles = (coords0+coords1+coords2)/3
                averageFiber = (p1Fiber+p2Fiber)/2

                # Find distance to nearest surface triangle
                distances = np.linalg.norm(averageTriangles-p2Fiber,axis=1)
                nearest = np.where(distances == np.amin(distances))

                # Store the plane of this triangle
                p0 = coords0[nearest,:]
                p1 = coords1[nearest,:]
                p2 = coords2[nearest,:]

                p01 = p1-p0
                p02 = p2-p0

                fiberVector = p2Fiber-p1Fiber

                # Compute distance to cutting plane
                t = (np.dot(np.squeeze(np.cross(p01,p02)),np.squeeze((p1Fiber-p0))))/(np.dot(np.squeeze(-fiberVector),np.squeeze(np.cross(p01,p02))))

                # New point 2 for fiber after cutting
                p2Fiber = p1Fiber+fiberVector*t

                # Obtain extents for floating bin
                binMin = np.amin(np.vstack((p1Fiber,p2Fiber)), axis=0)-maxAggD-np.linalg.norm(p1Fiber-p2Fiber)
                binMax = np.amax(np.vstack((p1Fiber,p2Fiber)), axis=0)+maxAggD+np.linalg.norm(p1Fiber-p2Fiber)

                # Verfiy cut fiber is inside the mesh      
                inside = False   
                inside = insideCheckFiber(vertices,tets,p1Fiber,p1Fiber+0.99999*fiberVector*t,\
                    binMin,binMax,coord1,coord2,coord3,coord4,maxC)

                if np.logical_and(inside == True,np.linalg.norm(p1Fiber-p2Fiber)<lFiber):

                    fiberLength = np.linalg.norm(p1Fiber-p2Fiber)
                    
                    return p1Fiber, p2Fiber, orienFiber, lFiber


            # If not trimming then discard fiber and try again
            else:

                pass

# Checks whether second fiber node falls within a tet (thus inside). 
def insideCheckFiber(vertices,tets,p1Fiber,p2Fiber,binMin,binMax,\
    coord1,coord2,coord3,coord4,maxC):

    # Store tet vertices that fall inside the bin
    coord1 = np.all([(coord1[:,0] > binMin[0]) , (coord1[:,0] < binMax[0]),\
        (coord1[:,1] > binMin[1]) , (coord1[:,1] < binMax[1]) ,\
        (coord1[:,2] > binMin[2]) , (coord1[:,2] < binMax[2])],axis=0)      
    coord2 = np.all([(coord2[:,0] > binMin[0]) , (coord2[:,0] < binMax[0]),\
        (coord2[:,1] > binMin[1]) , (coord2[:,1] < binMax[1]) ,\
        (coord2[:,2] > binMin[2]) , (coord2[:,2] < binMax[2])],axis=0)          
    coord3 = np.all([(coord3[:,0] > binMin[0]) , (coord3[:,0] < binMax[0]),\
        (coord3[:,1] > binMin[1]) , (coord3[:,1] < binMax[1]) ,\
        (coord3[:,2] > binMin[2]) , (coord3[:,2] < binMax[2])],axis=0)  
    coord4 = np.all([(coord4[:,0] > binMin[0]) , (coord4[:,0] < binMax[0]),\
        (coord4[:,1] > binMin[1]) , (coord4[:,1] < binMax[1]) ,\
        (coord4[:,2] > binMin[2]) , (coord4[:,2] < binMax[2])],axis=0)  

    binTets = np.any([coord1,coord2,coord3,coord4],axis=0)

    coord1 = vertices[tets.astype(int)[binTets,0]-1]
    coord2 = vertices[tets.astype(int)[binTets,1]-1]
    coord3 = vertices[tets.astype(int)[binTets,2]-1]
    coord4 = vertices[tets.astype(int)[binTets,3]-1]

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

    D10 = np.rot90(np.dstack((emptyOnes*p1Fiber[0],\
        emptyOnes*p1Fiber[1],emptyOnes*p1Fiber[2],emptyOnes)), 3)
    
    D1 = np.linalg.det(np.hstack((D10,D01,D02,D03)))
    D2 = np.linalg.det(np.hstack((D00,D10,D02,D03)))
    D3 = np.linalg.det(np.hstack((D00,D01,D10,D03)))
    D4 = np.linalg.det(np.hstack((D00,D01,D02,D10)))

    p1 = False

    if np.logical_and(np.logical_and(np.sign(D0) == np.sign(D1),\
        np.sign(D0) == np.sign(D2)),\
        np.logical_and(np.sign(D0) == np.sign(D3),\
        np.sign(D0) == np.sign(D4))).any():
        p1 = True

    D10 = np.rot90(np.dstack((emptyOnes*p2Fiber[0],\
        emptyOnes*p2Fiber[1],emptyOnes*p2Fiber[2],emptyOnes)), 3)
    
    D1 = np.linalg.det(np.hstack((D10,D01,D02,D03)))
    D2 = np.linalg.det(np.hstack((D00,D10,D02,D03)))
    D3 = np.linalg.det(np.hstack((D00,D01,D10,D03)))
    D4 = np.linalg.det(np.hstack((D00,D01,D02,D10)))

    p2 = False

    if np.logical_and(np.logical_and(np.sign(D0) == np.sign(D1),\
        np.sign(D0) == np.sign(D2)),\
        np.logical_and(np.sign(D0) == np.sign(D3),\
        np.sign(D0) == np.sign(D4))).any():
        p2 = True

    if np.logical_and(p1 == True,p2==True):
        return True