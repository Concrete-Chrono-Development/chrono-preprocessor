import numpy as np







def genTesselation(allNodes,allTets,aggDiameter,minAggD,geoName):
    
    # Create diameters list (including fictitious edge particle diameters)
    allDiameters = np.concatenate((np.array([1.1*minAggD,]*\
        int(len(allNodes)-len(aggDiameter))),aggDiameter))

    # Definition of Edge Points [Coordinates1,....,Coordinates6]
    edges1 = [allTets[:,0],allTets[:,1]]
    edges2 = [allTets[:,0],allTets[:,2]]
    edges3 = [allTets[:,0],allTets[:,3]]
    edges4 = [allTets[:,1],allTets[:,2]]
    edges5 = [allTets[:,1],allTets[:,3]]
    edges6 = [allTets[:,2],allTets[:,3]]

    edges = np.concatenate((edges1,edges2,edges3,edges4,edges5,edges6),\
        axis=1).T

    edges = np.sort(edges,axis=1)

    edgeNode1 = allNodes[(edges[:,0]-1).astype(int),:]
    edgeNode2 = allNodes[(edges[:,1]-1).astype(int),:]

    nodalDistance = np.linalg.norm(edgeNode1-edgeNode2, axis=1)
    edgeDistance = (nodalDistance - (allDiameters[(edges[:,0]-1).\
        astype(int)]/2) - (allDiameters[(edges[:,1]-1).astype(int)])/2)/2


    # Make unit vector from edgeNode 1 to edgeNode2, multiply vector 
    # by sum(agg1 and edgeDistance) and add to edgeNode2
    edgePoints = (edgeNode1-edgeNode2)/np.array([nodalDistance,]*3).T*\
        (np.array([allDiameters[(edges[:,1]-1).astype(int)],]*3).T/2+\
        np.array([edgeDistance,]*3).T)+edgeNode2


    # Form Edge Point List
    edgePoints = np.concatenate(np.split(edgePoints,6),axis=1)

    # Definition of Face Points faceNPoint[Coordinates]

    # Face 0 (Nodes: 1,2,3) 
    faceNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,15:18], axis=1),\
        np.linalg.norm(allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,12:15], axis=1),\
        np.linalg.norm(allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,9:12], axis=1))).T
    
    faceOffsetDistance = np.vstack(((faceNodalDistance[:,0] - allDiameters[(allTets[:,1]-1).astype(int)]/2)/2,\
        (faceNodalDistance[:,1] - allDiameters[(allTets[:,2]-1).astype(int)]/2)/2,\
        (faceNodalDistance[:,2] - allDiameters[(allTets[:,3]-1).astype(int)]/2)/2)).T

    facePoints = np.hstack(((allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,15:18])/\
        np.array([faceNodalDistance[:,0],]*3).T\
        *np.array([faceOffsetDistance[:,0],]*3).T + edgePoints[:,15:18],\
        (allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,12:15])/\
        np.array([faceNodalDistance[:,1],]*3).T\
        *np.array([faceOffsetDistance[:,1],]*3).T + edgePoints[:,12:15],\
        (allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,9:12])/\
        np.array([faceNodalDistance[:,2],]*3).T\
        *np.array([faceOffsetDistance[:,2],]*3).T + edgePoints[:,9:12]))

    face1x = (facePoints[:,0]+facePoints[:,3]+facePoints[:,6])/3
    face1y = (facePoints[:,1]+facePoints[:,4]+facePoints[:,7])/3
    face1z = (facePoints[:,2]+facePoints[:,5]+facePoints[:,8])/3

    face0Point = np.vstack((face1x,face1y,face1z)).T

    # Face 1 (Nodes: 0,2,3)
    faceNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,3:6], axis=1),\
        np.linalg.norm(allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,6:9], axis=1),\
        np.linalg.norm(allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,15:18], axis=1))).T
    
    faceOffsetDistance = np.vstack(((faceNodalDistance[:,0] - allDiameters[(allTets[:,3]-1).astype(int)]/2)/2,\
        (faceNodalDistance[:,1] - allDiameters[(allTets[:,2]-1).astype(int)]/2)/2,\
        (faceNodalDistance[:,2] - allDiameters[(allTets[:,0]-1).astype(int)]/2)/2)).T

    facePoints = np.hstack(((allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,3:6])/\
        np.array([faceNodalDistance[:,0],]*3).T\
        *np.array([faceOffsetDistance[:,0],]*3).T + edgePoints[:,3:6],\
        (allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,6:9])/\
        np.array([faceNodalDistance[:,1],]*3).T\
        *np.array([faceOffsetDistance[:,1],]*3).T + edgePoints[:,6:9],\
        (allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,15:18])/\
        np.array([faceNodalDistance[:,2],]*3).T\
        *np.array([faceOffsetDistance[:,2],]*3).T + edgePoints[:,15:18]))

    face1x = (facePoints[:,0]+facePoints[:,3]+facePoints[:,6])/3
    face1y = (facePoints[:,1]+facePoints[:,4]+facePoints[:,7])/3
    face1z = (facePoints[:,2]+facePoints[:,5]+facePoints[:,8])/3

    face1Point = np.vstack((face1x,face1y,face1z)).T

    # Face 2 (Nodes: 0,1,3)
    faceNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,12:15], axis=1),\
        np.linalg.norm(allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,6:9], axis=1),\
        np.linalg.norm(allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,0:3], axis=1))).T
    
    faceOffsetDistance = np.vstack(((faceNodalDistance[:,0] - allDiameters[(allTets[:,0]-1).astype(int)]/2)/2,\
        (faceNodalDistance[:,1] - allDiameters[(allTets[:,1]-1).astype(int)]/2)/2,\
        (faceNodalDistance[:,2] - allDiameters[(allTets[:,3]-1).astype(int)]/2)/2)).T

    facePoints = np.hstack(((allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,12:15])/\
        np.array([faceNodalDistance[:,0],]*3).T\
        *np.array([faceOffsetDistance[:,0],]*3).T + edgePoints[:,12:15],\
        (allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,6:9])/\
        np.array([faceNodalDistance[:,1],]*3).T\
        *np.array([faceOffsetDistance[:,1],]*3).T + edgePoints[:,6:9],\
        (allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,0:3])/\
        np.array([faceNodalDistance[:,2],]*3).T\
        *np.array([faceOffsetDistance[:,2],]*3).T + edgePoints[:,0:3]))

    face1x = (facePoints[:,0]+facePoints[:,3]+facePoints[:,6])/3
    face1y = (facePoints[:,1]+facePoints[:,4]+facePoints[:,7])/3
    face1z = (facePoints[:,2]+facePoints[:,5]+facePoints[:,8])/3

    face2Point = np.vstack((face1x,face1y,face1z)).T

    # Face 3 (Nodes: 0,1,2)
    faceNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,0:3], axis=1),\
        np.linalg.norm(allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,3:6], axis=1),\
        np.linalg.norm(allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,9:12], axis=1))).T
    
    faceOffsetDistance = np.vstack(((faceNodalDistance[:,0] - allDiameters[(allTets[:,2]-1).astype(int)]/2)/2,\
        (faceNodalDistance[:,1] - allDiameters[(allTets[:,1]-1).astype(int)]/2)/2,\
        (faceNodalDistance[:,2] - allDiameters[(allTets[:,0]-1).astype(int)]/2)/2)).T

    facePoints = np.hstack(((allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,0:3])/\
        np.array([faceNodalDistance[:,0],]*3).T\
        *np.array([faceOffsetDistance[:,0],]*3).T + edgePoints[:,0:3],\
        (allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,3:6])/\
        np.array([faceNodalDistance[:,1],]*3).T\
        *np.array([faceOffsetDistance[:,1],]*3).T + edgePoints[:,3:6],\
        (allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,9:12])/\
        np.array([faceNodalDistance[:,2],]*3).T\
        *np.array([faceOffsetDistance[:,2],]*3).T + edgePoints[:,9:12]))

    face1x = (facePoints[:,0]+facePoints[:,3]+facePoints[:,6])/3
    face1y = (facePoints[:,1]+facePoints[:,4]+facePoints[:,7])/3
    face1z = (facePoints[:,2]+facePoints[:,5]+facePoints[:,8])/3

    face3Point = np.vstack((face1x,face1y,face1z)).T

    # Definition of Tet-Points [Coordinates]
    tetNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,0]-1).astype(int),:]-face0Point, axis=1),\
        np.linalg.norm(allNodes[(allTets[:,1]-1).astype(int),:]-face1Point, axis=1),\
        np.linalg.norm(allNodes[(allTets[:,2]-1).astype(int),:]-face2Point, axis=1),\
        np.linalg.norm(allNodes[(allTets[:,3]-1).astype(int),:]-face3Point, axis=1))).T     

    tetOffsetDistance = np.vstack(((tetNodalDistance[:,0] - allDiameters[(allTets[:,0]-1).astype(int)]/2)/2,\
        (tetNodalDistance[:,1] - allDiameters[(allTets[:,1]-1).astype(int)]/2)/2,\
        (tetNodalDistance[:,2] - allDiameters[(allTets[:,2]-1).astype(int)]/2)/2,\
        (tetNodalDistance[:,3] - allDiameters[(allTets[:,3]-1).astype(int)]/2)/2)).T

    tetPoints = np.hstack(((allNodes[(allTets[:,0]-1).astype(int),:]-face0Point)/\
        np.array([tetNodalDistance[:,0],]*3).T*\
        np.array([tetOffsetDistance[:,0],]*3).T + face0Point,\
        (allNodes[(allTets[:,1]-1).astype(int),:]-face1Point)/\
        np.array([tetNodalDistance[:,1],]*3).T*\
        np.array([tetOffsetDistance[:,1],]*3).T + face1Point,\
        (allNodes[(allTets[:,2]-1).astype(int),:]-face2Point)/\
        np.array([tetNodalDistance[:,2],]*3).T*\
        np.array([tetOffsetDistance[:,2],]*3).T + face2Point,\
        (allNodes[(allTets[:,3]-1).astype(int),:]-face3Point)/\
        np.array([tetNodalDistance[:,3],]*3).T*\
        np.array([tetOffsetDistance[:,3],]*3).T + face3Point))

    tetx = (tetPoints[:,0]+tetPoints[:,3]+tetPoints[:,6]+tetPoints[:,9])/4
    tety = (tetPoints[:,1]+tetPoints[:,4]+tetPoints[:,7]+tetPoints[:,10])/4
    tetz = (tetPoints[:,2]+tetPoints[:,5]+tetPoints[:,8]+tetPoints[:,11])/4

    tetPoints = np.vstack((tetx,tety,tetz)).T

    # Coordinates for each facet in each tet
    facet1 = np.concatenate(([tetPoints,face0Point,edgePoints[:,9:12]]),axis=1)
    facet2 = np.concatenate(([tetPoints,face0Point,edgePoints[:,12:15]]),axis=1)
    facet3 = np.concatenate(([tetPoints,face0Point,edgePoints[:,15:18]]),axis=1)
    facet4 = np.concatenate(([tetPoints,face1Point,edgePoints[:,3:6]]),axis=1)
    facet5 = np.concatenate(([tetPoints,face1Point,edgePoints[:,6:9]]),axis=1)
    facet6 = np.concatenate(([tetPoints,face1Point,edgePoints[:,15:18]]),axis=1)
    facet7 = np.concatenate(([tetPoints,face2Point,edgePoints[:,0:3]]),axis=1)
    facet8 = np.concatenate(([tetPoints,face2Point,edgePoints[:,6:9]]),axis=1)
    facet9 = np.concatenate(([tetPoints,face2Point,edgePoints[:,12:15]]),axis=1)
    facet10 = np.concatenate(([tetPoints,face3Point,edgePoints[:,0:3]]),axis=1)
    facet11 = np.concatenate(([tetPoints,face3Point,edgePoints[:,3:6]]),axis=1)
    facet12 = np.concatenate(([tetPoints,face3Point,edgePoints[:,9:12]]),axis=1)

    # Facets for P0 in General Tet [facet[Point1[x,y,z],Point2[x,y,z],Point3[x,y,z]],....,facetN]
    facetsP0 = np.concatenate(([facet4,facet5,facet7,facet8,facet10,facet11]),axis=1)

    # Facets for P1 in General Tet [facet[Point1[x,y,z],Point2[x,y,z],Point3[x,y,z]],....,facetN]
    facetsP1 = np.concatenate(([facet1,facet2,facet7,facet9,facet10,facet12]),axis=1)

    # Facets for P2 in General Tet [facet[Point1[x,y,z],Point2[x,y,z],Point3[x,y,z]],....,facetN]
    facetsP2 = np.concatenate(([facet1,facet3,facet4,facet6,facet11,facet12]),axis=1)

    # Facets for P3 in General Tet [facet[Point1[x,y,z],Point2[x,y,z],Point3[x,y,z]],....,facetN]
    facetsP3 = np.concatenate(([facet2,facet3,facet5,facet6,facet8,facet9]),axis=1)

    # Combination of all facets for general tet (24)
    cellTetFacets = np.concatenate(([facetsP0,facetsP1,facetsP2,facetsP3]),\
        axis=1)

    # Combination of nonrepeating facets for general tet (12)
    tetFacets = np.concatenate(([facet1,facet2,facet3,facet4,facet5,facet6,\
        facet7,facet8,facet9,facet10,facet11,facet12]),axis=1)
    facets = tetFacets.reshape(-1,9)

    # Initialize cell facet list
    #cells = []

    #for i in range(1,len(allNodes)):

        # Find locations in tet list which contain the given vertex
    #    location = np.argwhere(allTets==i)

        # Initialize empty facet list for vertex
    #    cellFacets=np.zeros(([len(location)*6,9]))

        # Store facets for given cell
    #    for x in range(0,len(location)):
    #       cellFacets[x*6:x*6+6,:] = cellTetFacets[location[x,0],\
    #        location[x,1]*54:location[x,1]*54+54].reshape((-1,9))

        # Add cell facets to master list
    #    cells.append(cellFacets)

    threes = 3*np.array([np.ones((len(facets)))]).T

    facetCenters = np.array((facets[:,0]+facets[:,3]+facets[:,6],\
        facets[:,1]+facets[:,4]+facets[:,7],\
        facets[:,2]+facets[:,5]+facets[:,8])).T/threes

    vectorAB = np.array((tetFacets.reshape(-1,9)[:,3]-\
        tetFacets.reshape(-1,9)[:,0],tetFacets.reshape(-1,9)[:,4]-\
        tetFacets.reshape(-1,9)[:,1],tetFacets.reshape(-1,9)[:,5]-\
        tetFacets.reshape(-1,9)[:,2])).T
    vectorAC = np.array((tetFacets.reshape(-1,9)[:,6]-\
        tetFacets.reshape(-1,9)[:,0],tetFacets.reshape(-1,9)[:,7]-\
        tetFacets.reshape(-1,9)[:,1],tetFacets.reshape(-1,9)[:,8]-\
        tetFacets.reshape(-1,9)[:,2])).T

    facetAreas = 0.5 * np.linalg.norm(np.cross(vectorAB,vectorAC),axis=1)
    facetNormals = np.cross(vectorAB,vectorAC)/np.array([np.linalg.norm\
        (np.cross(vectorAB,vectorAC),axis=1),]*3).T

    # Specify tet-node connectivity for facets (i.e. facet 1 connected by 
        # node 1 and 2)
    tetn1 = [1,1,2,0,0,2,0,0,1,0,0,1]
    tetn2 = [2,3,3,2,3,3,1,3,3,1,2,2]


    return tetFacets, facetCenters, facetAreas, facetNormals, \
        tetn1, tetn2, tetPoints, allDiameters





