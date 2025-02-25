## ===========================================================================
## CHRONO WORKBENCH:github.com/Concrete-Chrono-Development/chrono-preprocessor
##
## Copyright (c) 2024 
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
## Calculate the intersection of fiber with facets
##
## ===========================================================================


import numpy as np
import math

def gen_LDPMCSL_facetfiberInt(p1Fiber,p2Fiber,dFiber,lFiber,orienFibers,\
        geoName,allTets,allNodes,tetFacets,dataList,tetn1,tetn2,facetNormals,facetCenters):
    
    # Number of total fiber
    NumberofFibers = len(p1Fiber[:,1])
    #print('NumberofFibers',NumberofFibers)
    # Initialize a data array and matrix for interseted fiber
    FiberdataList = []
    ProjectedFacet = []
    FibertetList = []
    IntersectedFiber = []
    No = 0.0
    TotalTet = 0.0
    TotalFiber = 0.0


    for z in range(0,NumberofFibers):
        # Obtain extents for floating bin for fiber
        binMin = np.amin(np.vstack((p1Fiber[z,:],p2Fiber[z,:])), axis=0)-1.0*lFiber[z]-dFiber
        binMax = np.amax(np.vstack((p1Fiber[z,:],p2Fiber[z,:])), axis=0)+1.0*lFiber[z]+dFiber

        for x in range(0,len(allTets)):
            # Store tet vertices that fall inside the bin

            coord1 = allNodes[(allTets[x,0]-1).astype(int)]
            coord2 = allNodes[(allTets[x,1]-1).astype(int)]
            coord3 = allNodes[(allTets[x,2]-1).astype(int)]
            coord4 = allNodes[(allTets[x,3]-1).astype(int)]

            coord1 = np.all([(coord1[0] > binMin[0]) , (coord1[0] < binMax[0]),\
                (coord1[1] > binMin[1]) , (coord1[1] < binMax[1]) ,\
                (coord1[2] > binMin[2]) , (coord1[2] < binMax[2])],axis=0)      
            coord2 = np.all([(coord2[0] > binMin[0]) , (coord2[0] < binMax[0]),\
                (coord2[1] > binMin[1]) , (coord2[1] < binMax[1]) ,\
                (coord2[2] > binMin[2]) , (coord2[2] < binMax[2])],axis=0)          
            coord3 = np.all([(coord3[0] > binMin[0]) , (coord3[0] < binMax[0]),\
                (coord3[1] > binMin[1]) , (coord3[1] < binMax[1]) ,\
                (coord3[2] > binMin[2]) , (coord3[2] < binMax[2])],axis=0)  
            coord4 = np.all([(coord4[0] > binMin[0]) , (coord4[0] < binMax[0]),\
                (coord4[1] > binMin[1]) , (coord4[1] < binMax[1]) ,\
                (coord4[2] > binMin[2]) , (coord4[2] < binMax[2])],axis=0)  

            binTets = np.any([coord1,coord2,coord3,coord4],axis=0) 
            
            if binTets == True:
                FibertetList.append(np.array([int(x),int(z)]))  

    FibertetList=np.array(FibertetList).reshape(-1,2)
    FibertetList=FibertetList[FibertetList[:,0].argsort()]                          
    FiberBin=np.unique(FibertetList[:,1])
    #print('FiberBin',len(FiberBin))

    #######################################################
    # For fiberfacet using Matthew's rotational matrix

    # Store facets as set of three points
    facets = tetFacets.reshape(-1,9)

    # Store particles accociated with facets
    p1 = allNodes[(allTets[:,tetn1]-1).astype(int),:].reshape(-1,3)
    p2 = allNodes[(allTets[:,tetn2]-1).astype(int),:].reshape(-1,3)


    # Projected facet normal
    pn = (p2-p1)/np.array([np.linalg.norm(p2-p1,axis=1),]*3).T
    pn = pn.reshape(-1,3)
    #print('n',pn)

    InitialNormal=[]
    coords=[]

    for x in range(0,len(allTets)):

        for y in range(0,12):


            Check = pn[12*x+y,0]*facetNormals[12*x+y,0]+\
                    pn[12*x+y,1]*facetNormals[12*x+y,1]+\
                    pn[12*x+y,2]*facetNormals[12*x+y,2]

            if Check < 0.0:
                InitialN = -1.0*facetNormals[12*x+y,0:3]
                c1 = facets[12*x+y,0:3]
                c2 = facets[12*x+y,6:9]
                c3 = facets[12*x+y,3:6]

            else:
                InitialN = facetNormals[12*x+y,0:3]
                c1 = facets[12*x+y,0:3]
                c2 = facets[12*x+y,3:6]
                c3 = facets[12*x+y,6:9]

            coords.append(np.array([c1,c2,c3]))
            InitialNormal.append(InitialN)

    InitialNormal = np.array(InitialNormal).reshape(-1,3)
    coords = np.array(coords).reshape(-1,9)
    facetNormals = InitialNormal
    facets = coords.reshape(-1,9)

    # Formation of rotation stacked matrix (3 x 3 x nFacets)
    v = np.cross(facetNormals,pn.reshape(-1,3))
    zeros = np.zeros(len(v),)
    ssc = np.array(([[zeros, -v[:,2], v[:,1]],[ v[:,2], zeros, -v[:,0]],\
        [ -v[:,1], v[:,0], zeros]]))
    identity = np.dstack([np.eye(3)]*len(v))
    mulNormalsPn = np.matmul(np.expand_dims(facetNormals, axis=2),\
        np.expand_dims(pn.reshape(-1,3), axis=1)).T
    R = identity + ssc + (np.matmul(ssc.T,ssc.T).T)*(1-mulNormalsPn)/\
    (np.dot(np.linalg.norm(v),np.linalg.norm(v)))

    # Make vectors from center to corners of the facets
    vectorOA = np.expand_dims(facetCenters[:,0:3]-facets[:,0:3], axis=2)
    vectorOB = np.expand_dims(facetCenters[:,0:3]-facets[:,3:6], axis=2)
    vectorOC = np.expand_dims(facetCenters[:,0:3]-facets[:,6:9], axis=2)

    pvectorOA =np.squeeze(np.matmul(np.transpose(R.T,(0, 2, 1)),vectorOA))
    pvectorOB = np.squeeze(np.matmul(np.transpose(R.T,(0, 2, 1)),vectorOB))
    pvectorOC = np.squeeze(np.matmul(np.transpose(R.T,(0, 2, 1)),vectorOC))

    RcoordP1= facetCenters-pvectorOA
    RcoordP2= facetCenters-pvectorOB
    RcoordP3= facetCenters-pvectorOC

    # Compute the plane supporting the triangle (coordP1, coordP2, coordP3)  normal: nprojectedfacet offset: d
    vector12 = RcoordP2-RcoordP1
    vector23 = RcoordP3-RcoordP2
    vector13 = RcoordP3-RcoordP1
    vector31 = RcoordP1-RcoordP3

    Normal = np.cross(vector12,vector13)/np.array([np.linalg.norm(np.cross(vector12,vector13),axis=1),]*3).T

    RcoordP1=np.array(RcoordP1).reshape(-1,3)
    RcoordP2=np.array(RcoordP2).reshape(-1,3)
    RcoordP3=np.array(RcoordP3).reshape(-1,3)
    vector12=np.array(vector12).reshape(-1,3)
    vector31=np.array(vector31).reshape(-1,3)
    vector23=np.array(vector23).reshape(-1,3)
    vector13=np.array(vector13).reshape(-1,3)
    Normal=np.array(Normal).reshape(-1,3)
    projectedFacet = np.concatenate((RcoordP1,RcoordP2,RcoordP3))

    for i in range(0,len(FibertetList)):

        x=int(FibertetList[i,0])
        z=int(FibertetList[i,1])

        TotalTet = TotalTet + 1    

        # Check the intersection of this inside fiber with 12 facsets of tet
        for y in range(0,12):

            p12Fiber = np.array([p2Fiber[z,0]-p1Fiber[z,0],\
                p2Fiber[z,1]-p1Fiber[z,1],\
                p2Fiber[z,2]-p1Fiber[z,2]])


            offsetd = Normal[12*x+y,0]*RcoordP1[12*x+y,0]+\
                Normal[12*x+y,1]*RcoordP1[12*x+y,1]+\
                Normal[12*x+y,2]*RcoordP1[12*x+y,2]

            offsetd = -1*offsetd

            facetnormalDotp12Fiber = Normal[12*x+y,0]*p12Fiber[0]+\
                Normal[12*x+y,1]*p12Fiber[1]+\
                Normal[12*x+y,2]*p12Fiber[2]
            

            facetnormalDotp1Fiber = Normal[12*x+y,0]*p1Fiber[z,0]+\
                Normal[12*x+y,1]*p1Fiber[z,1]+\
                Normal[12*x+y,2]*p1Fiber[z,2]

            # Ignore line parallel to (or lying in) the plane
            if abs(facetnormalDotp12Fiber) > 0.0:
                t = - (offsetd + facetnormalDotp1Fiber)/facetnormalDotp12Fiber

            # Check if the intersection point is between p1Fiber and p2Fiber
                if t >= 0.0 and t <= 1.0:

                    # From here can find why we have loop over tets, facets and fibers
                    P = np.array([p1Fiber[z,0] + t*p12Fiber[0],\
                        p1Fiber[z,1] + t*p12Fiber[1],\
                        p1Fiber[z,2] + t*p12Fiber[2]])

                    allPoints =[]
                    CheckEdge1 = []
                    CheckEdge2 = []
                    CheckEdge3 = []

                    # N12 = np.cross(Normal[12*x+y,0:3],vector12[12*x+y,0:3])
                    # PP1 = np.array([P[0]-RcoordP1[12*x+y,0],P[1]-RcoordP1[12*x+y,1],P[2]-RcoordP1[12*x+y,2]])
                    # CheckEdge1 = np.dot(PP1.T,N12)/np.linalg.norm(N12)

                    PP1 = np.array([P[0]-RcoordP1[12*x+y,0],P[1]-RcoordP1[12*x+y,1],P[2]-RcoordP1[12*x+y,2]])
                    N12 = np.cross(vector12[12*x+y,0:3],PP1)
                    CheckEdge1 = np.dot(Normal[12*x+y,0:3].T,N12)

                    if (CheckEdge1>=0).any():

                        # N23 = np.cross(Normal[12*x+y,0:3],vector23[12*x+y,0:3])
                        # PP2 = np.array([P[0]-RcoordP2[12*x+y,0],P[1]-RcoordP2[12*x+y,1],P[2]-RcoordP2[12*x+y,2]])
                        # CheckEdge2 = np.dot(PP2.T,N23)/np.linalg.norm(N23)

                        PP2 = np.array([P[0]-RcoordP2[12*x+y,0],P[1]-RcoordP2[12*x+y,1],P[2]-RcoordP2[12*x+y,2]])
                        N23 = np.cross(vector23[12*x+y,0:3],PP2)
                        CheckEdge2 = np.dot(Normal[12*x+y,0:3].T,N23)

                        if (CheckEdge2>=0).any():

                            # N31 = np.cross(Normal[12*x+y,0:3],vector31[12*x+y,0:3])
                            # PP3 = np.array([P[0]-RcoordP3[12*x+y,0],P[1]-RcoordP3[12*x+y,1],P[2]-RcoordP3[12*x+y,2]])
                            # CheckEdge3 = np.dot(PP3.T,N31)/np.linalg.norm(N31)

                            PP3 = np.array([P[0]-RcoordP3[12*x+y,0],P[1]-RcoordP3[12*x+y,1],P[2]-RcoordP3[12*x+y,2]])
                            N31 = np.cross(vector31[12*x+y,0:3],PP3)
                            CheckEdge3 = np.dot(Normal[12*x+y,0:3].T,N31)

                            if (CheckEdge3>=0).any():

                                # Determination of Short and Long lenght of fiber intersected facet
                                IntersectedFiber.append(np.array([int(z)]))
                                No = No + 1
                                distancep1Fiber = math.sqrt(((P[0]-p1Fiber[z,0])**2)+\
                                    ((P[1]-p1Fiber[z,1])**2)+\
                                    ((P[2]-p1Fiber[z,2])**2))
                                distancep2Fiber = math.sqrt(((P[0]-p2Fiber[z,0])**2)+\
                                    ((P[1]-p2Fiber[z,1])**2)+\
                                    ((P[2]-p2Fiber[z,2])**2))

                                tetindex = x+1
                                facetindex = y+1
                                FiberShortLenght = distancep2Fiber
                                FiberLongLenght = distancep1Fiber

                                if  distancep1Fiber < distancep2Fiber:
                                    FiberShortLenght = distancep1Fiber
                                    FiberLongLenght = distancep2Fiber

                                InterPerFacet = 1.0

                                OneIntersectedFiber = np.array([tetindex, facetindex, orienFibers[z,0], orienFibers[z,1], orienFibers[z,2],\
                                    FiberShortLenght, FiberLongLenght, dFiber, InterPerFacet])

                                FiberdataList.append(OneIntersectedFiber)
            TotalIntersections = No

    FiberdataList = np.array(FiberdataList)
    FiberdataList = FiberdataList.reshape(-1,9)
    FiberdataList=FiberdataList[FiberdataList[:,0].argsort()]
    FiberdataList=FiberdataList[np.lexsort(FiberdataList[:,::-1].T)]
    IntersectedFiber = np.unique(np.array(IntersectedFiber))
    TotalFiber = len(IntersectedFiber)
    #print('TotalFiber',TotalFiber)
    #print(IntersectedFiber)


    for k in range(1,len(FiberdataList)):
        if ((FiberdataList[k,0]==FiberdataList[k-1,0]) and (FiberdataList[k,1]==FiberdataList[k-1,1])):
            FiberdataList[k,8] = FiberdataList[k-1,8] + 1

    MaxInterPerFacet = np.max(FiberdataList[:,8])

    return FiberdataList,TotalIntersections,MaxInterPerFacet,TotalTet,TotalFiber,IntersectedFiber,projectedFacet    