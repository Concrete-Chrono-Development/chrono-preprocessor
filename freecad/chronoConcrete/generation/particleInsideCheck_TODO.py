import numpy as np

def insideCheck(vertices,tets,center,parDiameter,binMin,binMax,coord1,\
    coord2,coord3,coord4,maxC):

    # Store tet vertices that fall inside the bin
    binTets = np.any([(coord1[:,0] > binMin[0]) & (coord1[:,0] < binMax[0]) &
                      (coord1[:,1] > binMin[1]) & (coord1[:,1] < binMax[1]) &
                      (coord1[:,2] > binMin[2]) & (coord1[:,2] < binMax[2]),
                      (coord2[:,0] > binMin[0]) & (coord2[:,0] < binMax[0]) &
                      (coord2[:,1] > binMin[1]) & (coord2[:,1] < binMax[1]) &
                      (coord2[:,2] > binMin[2]) & (coord2[:,2] < binMax[2]),
                      (coord3[:,0] > binMin[0]) & (coord3[:,0] < binMax[0]) &
                      (coord3[:,1] > binMin[1]) & (coord3[:,1] < binMax[1]) &
                      (coord3[:,2] > binMin[2]) & (coord3[:,2] < binMax[2]),
                      (coord4[:,0] > binMin[0]) & (coord4[:,0] < binMax[0]) &
                      (coord4[:,1] > binMin[1]) & (coord4[:,1] < binMax[1]) &
                      (coord4[:,2] > binMin[2]) & (coord4[:,2] < binMax[2])], axis=0)


    node = np.array([    [center[0,0]+parDiameter/2, center[0,1]+parDiameter/2, center[0,2]+parDiameter/2],
        [center[0,0]-parDiameter/2, center[0,1]+parDiameter/2, center[0,2]+parDiameter/2],
        [center[0,0]+parDiameter/2, center[0,1]-parDiameter/2, center[0,2]+parDiameter/2],
        [center[0,0]+parDiameter/2, center[0,1]+parDiameter/2, center[0,2]-parDiameter/2],
        [center[0,0]-parDiameter/2, center[0,1]-parDiameter/2, center[0,2]+parDiameter/2],
        [center[0,0]+parDiameter/2, center[0,1]-parDiameter/2, center[0,2]-parDiameter/2],
        [center[0,0]-parDiameter/2, center[0,1]+parDiameter/2, center[0,2]-parDiameter/2],
        [center[0,0]-parDiameter/2, center[0,1]-parDiameter/2, center[0,2]-parDiameter/2]
    ], dtype=np.float32)


    coord1 = vertices[tets.astype(int)[binTets,0]-1]
    coord2 = vertices[tets.astype(int)[binTets,1]-1]
    coord3 = vertices[tets.astype(int)[binTets,2]-1]
    coord4 = vertices[tets.astype(int)[binTets,3]-1]

    emptyOnes = np.ones(len(coord1[:,0]))

    D00 = np.rot90(np.dstack((coord1[:,0],coord1[:,1],coord1[:,2],emptyOnes)), 3)
    D01 = np.rot90(np.dstack((coord2[:,0],coord2[:,1],coord2[:,2],emptyOnes)), 3)
    D02 = np.rot90(np.dstack((coord3[:,0],coord3[:,1],coord3[:,2],emptyOnes)), 3)
    D03 = np.rot90(np.dstack((coord4[:,0],coord4[:,1],coord4[:,2],emptyOnes)), 3)

    D0 = np.linalg.det(np.concatenate((D00,D01,D02,D03), axis=1))

    D10 = np.rot90(np.dstack((emptyOnes * node[:,0], emptyOnes * node[:,1],
                             emptyOnes * node[:,2], emptyOnes)), 3)
    D1 = np.linalg.det(np.concatenate((D10,D01,D02,D03), axis=2))
    D2 = np.linalg.det(np.concatenate((D00,D10,D02,D03), axis=2))
    D3 = np.linalg.det(np.concatenate((D00,D01,D10,D03), axis=2))
    D4 = np.linalg.det(np.concatenate((D00,D01,D02,D10), axis=2))

    if np.count_nonzero(np.all(np.sign(D0) == np.sign(D1), axis=0) &
                            np.all(np.sign(D0) == np.sign(D2), axis=0) &
                            np.all(np.sign(D0) == np.sign(D3), axis=0) &
                            np.all(np.sign(D0) == np.sign(D4), axis=0)) == 8:
        return True
    else:
        return False