import numpy as np



def overlapCheck(vertices,center,aggDiameter,facePoints,binMin,binMax,\
    minPar,maxEdgeLength,aggOffset,parDiameterList):

    # Store particle vertices that fall inside the bin
    binTestParticles = np.all([(vertices[:,0] > binMin[0]) , \
        (vertices[:,0] < binMax[0]) , (vertices[:,1] > binMin[1]) , \
        (vertices[:,1] < binMax[1]) , (vertices[:,2] > binMin[2]) , \
        (vertices[:,2] < binMax[2])],axis=0)
    existingvertices = vertices[binTestParticles,:]
    existingAggD = parDiameterList[binTestParticles]

    # Compute distance between particles 
    if len(existingvertices>0):
        nodalDistance = np.linalg.norm(center-existingvertices, axis=1)
        aggOffsetDist = nodalDistance - aggDiameter/2 - existingAggD\
            /2 - aggOffset
    else: 
        aggOffsetDist = np.array(([1]))

    # Kill and return if overlap
    if (aggOffsetDist<0).any():
        return True,"NA"

    # Store edge vertices that fall inside the bin
    binTestSurf = np.all([(facePoints[:,0] > binMin[0]) , \
        (facePoints[:,0] < binMax[0]) , (facePoints[:,1] > \
            binMin[1]) , (facePoints[:,1] < binMax[1]) ,\
        (facePoints[:,2] > binMin[2]) , (facePoints[:,2] < \
            binMax[2])],axis=0)     
    existingSurf = facePoints[binTestSurf,:]

    # Compute distance between particle and edge vertices
    if len(existingSurf>0):
        surfNodalDistance = np.linalg.norm(center-existingSurf, axis=1)
        aggSurfaceDist = surfNodalDistance - aggDiameter/2 - 1.1*minPar/2
    else: 
        aggSurfaceDist = np.array(([1]))

    # Kill and return if overlap
    if (aggSurfaceDist<0).any():
        return True,"NA"

    # Otherwise return false and check if critically close to surface
    if len(existingSurf>0):
        if (surfNodalDistance<math.sqrt(1/3*maxEdgeLength**2+(aggDiameter/2)**2)).any():
            return False,True

    return False,False