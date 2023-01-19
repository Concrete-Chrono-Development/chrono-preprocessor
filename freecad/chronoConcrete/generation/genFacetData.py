import numpy as np








def genFacetData(allNodes,allTets,tetFacets,facetCenters,\
    facetAreas,facetNormals,tetn1,tetn2,materialList,materialRule,\
    multiMaterial,cementStructure,edgeMaterialList):

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


    InitialNormal=np.array(InitialNormal).reshape(-1,3)
    coords = np.array(coords).reshape(-1,9)
    facetNormals=InitialNormal
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

    # Clear not needed variables from memory
    del Check
    del c1
    del c2
    del c3
    del v
    del zeros
    del identity
    del mulNormalsPn
    del ssc        

    # Define 1st projected tangential
    tan1 = np.expand_dims(facetCenters-facets[:,0:3], axis=2)
    ptan1 = np.squeeze(np.matmul(np.transpose(R.T,(0, 2, 1)),tan1))/\
        np.array([np.linalg.norm(np.squeeze(np.matmul(np.transpose(R.T,\
            (0, 2, 1)),tan1)),axis=1),]*3).T

    # Define 2nd projected tangential
    ptan2 = np.cross(pn,ptan1)/np.array([np.linalg.norm(np.cross(pn,ptan1),\
        axis=1),]*3).T

    # Sub-tet Volume
    coord1 = facets[:,0:3]
    coord2 = facets[:,3:6]
    coord3 = facets[:,6:9]
    volCalc1 = np.expand_dims(np.transpose(coord1-p1).T, axis=1)
    volCalc2 = np.expand_dims(np.transpose(np.cross((coord2-p1),\
        (coord3-p1))).T, axis=2)
    volCalc3 = np.expand_dims(np.transpose(coord1-p2).T, axis=1)
    volCalc4 = np.expand_dims(np.transpose(np.cross((coord2-p2),\
        (coord3-p2))).T, axis=2)
    facetVol1 = np.squeeze(abs(np.matmul(volCalc1,volCalc2))/6)
    facetVol2 = np.squeeze(abs(np.matmul(volCalc3,volCalc4))/6)
    subtetVol = np.squeeze(abs(np.matmul(volCalc1,volCalc2))/6 \
        + abs(np.matmul(volCalc3,volCalc4))/6)
    #print(facetVol1)
    #print(facetVol2)

    # Clear not needed variables from memory
    del R
    del tan1
    del coord1
    del coord2
    del coord3
    del p1
    del p2
    del volCalc1
    del volCalc2
    del volCalc3
    del volCalc4

    # Projected Area
    areaCalc1 = np.squeeze(np.matmul(np.expand_dims(facetNormals, axis=1),\
        np.expand_dims(pn, axis=2)))
    areaCalc2 = np.linalg.norm(facetNormals, axis=1)\
        *np.linalg.norm(pn, axis=1)
    pArea = abs(areaCalc1/areaCalc2*facetAreas)

    # Initialize a data matrix for all facet data
    dataList = np.empty([len(allTets)*36,10])

    # Initialize a matrix for facet material information
    facetMaterial = np.empty([len(allTets)*12,])

    # Initialize a matrix for particle material information
    particleMaterial = np.empty([len(allTets)*12,2])


    # Extend material list for edge nodes
    if multiMaterial in ['on','On','Y','y','Yes','yes']:

        materialList = np.concatenate((2*np.ones([len(allNodes)-\
            len(materialList),]),materialList))

    elif cementStructure in ['on','On','Y','y','Yes','yes']:

    # materialList = np.concatenate((0*np.ones([len(allNodes)-\
    #     len(materialList),]),materialList))  
        materialList = np.concatenate((edgeMaterialList,materialList))

    else:

        materialList = np.concatenate((0*np.ones([len(allNodes)-\
            len(materialList),]),materialList))            

    for x in range(0,len(allTets)):

        for y in range(0,12):

            # Store: [n1 n2 cx cy cz fa nx ny nz vol]
            dataList[36*x+3*y,0] = allTets[x,tetn1[y]].astype(int)
            dataList[36*x+3*y,1] = allTets[x,tetn2[y]].astype(int)
            dataList[36*x+3*y,2:5] = facetCenters[12*x+y,:]
            dataList[36*x+3*y,5] = facetAreas[12*x+y]
            dataList[36*x+3*y,6:9] = facetNormals[12*x+y,:]
            dataList[36*x+3*y,9] = subtetVol[12*x+y]

            # Store: [pa px py pz qx qy qz sx sy sz]
            dataList[36*x+3*y+1,0] = pArea[12*x+y]
            dataList[36*x+3*y+1,1:4] = pn[12*x+y,:]
            dataList[36*x+3*y+1,4:7] = ptan1[12*x+y,:]
            dataList[36*x+3*y+1,7:10] = ptan2[12*x+y,:]

            # Store: [cx1 cy1 cz1 cx2 cy2 cz2 cx3 cy3 cz3 0]
        # if multiMaterial in ['on','On','Y','y','Yes','yes']
            # 1 = ITZ
            # 2 = Binder
            # 3 = Aggregate
            dataList[36*x+3*y+2,0:9] = facets[12*x+y,:]

            # Store particle materials
            particleMaterial[12*x+y,0] = materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)
            particleMaterial[12*x+y,1] = materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)

            # Material rule based on particle diameters volumes
            if materialRule == 9 or materialRule == 10:

                if facetVol1.all() == facetVol2.all():
                    
                    dataList[36*x+3*y+2,9] = materialList[allTets[x,tetn1[y]]\
                        .astype(int)-1].astype(int)
                    facetMaterial[12*x+y] = materialList[allTets[x,tetn1[y]]\
                        .astype(int)-1].astype(int)

                else:

                    if facetVol1.any() > facetVol2.any():

                        dataList[36*x+3*y+2,9] = materialList[allTets[x,tetn1[y]]\
                            .astype(int)-1].astype(int)
                        facetMaterial[12*x+y] = materialList[allTets[x,tetn1[y]]\
                            .astype(int)-1].astype(int)

                    else: 

                        dataList[36*x+3*y+2,9] = materialList[allTets[x,tetn2[y]]\
                            .astype(int)-1].astype(int)
                        facetMaterial[12*x+y] = materialList[allTets[x,tetn2[y]]\
                            .astype(int)-1].astype(int)

            elif materialRule > 0:

                if materialRule == 1:
                    aggITZ = 3
                    aggBinder = 3
                    itzBinder = 1
                elif materialRule == 2:
                    aggITZ = 3
                    aggBinder = 3
                    itzBinder = 2
                elif materialRule == 3:
                    aggITZ = 3
                    aggBinder = 2
                    itzBinder = 1
                elif materialRule == 4:
                    aggITZ = 3
                    aggBinder = 2
                    itzBinder = 2
                elif materialRule == 5:
                    aggITZ = 1
                    aggBinder = 3
                    itzBinder = 1
                elif materialRule == 6:
                    aggITZ = 1
                    aggBinder = 3
                    itzBinder = 2
                elif materialRule == 7:
                    aggITZ = 1
                    aggBinder = 2
                    itzBinder = 1
                elif materialRule == 8:
                    aggITZ = 1
                    aggBinder = 2
                    itzBinder = 2
                
                if materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int) \
                    == materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int):
                    
                    dataList[36*x+3*y+2,9] = materialList[allTets[x,tetn1[y]]\
                        .astype(int)-1].astype(int)
                    facetMaterial[12*x+y] = materialList[allTets[x,tetn1[y]]\
                        .astype(int)-1].astype(int)
                
                elif materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)==1 and materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)==2:

                    dataList[36*x+3*y+2,9] = itzBinder
                    facetMaterial[12*x+y] = itzBinder

                elif materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)==2 and materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)==1:

                    dataList[36*x+3*y+2,9] = itzBinder
                    facetMaterial[12*x+y] = itzBinder

                elif materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)==1 and materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)==3:

                    dataList[36*x+3*y+2,9] = aggITZ
                    facetMaterial[12*x+y] = aggITZ

                elif materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)==3 and materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)==1:

                    dataList[36*x+3*y+2,9] = aggITZ
                    facetMaterial[12*x+y] = aggITZ

                else:
                    
                    dataList[36*x+3*y+2,9] = aggBinder
                    facetMaterial[12*x+y] = aggBinder

            else:

                dataList[36*x+3*y+2,9] = 0
                facetMaterial[12*x+y] = 0                       

    return dataList,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial
