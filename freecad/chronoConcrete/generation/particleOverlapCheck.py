import math
import numpy as np



def overlapCheck(nodes,center,aggDiameter,facePoints,binMin,binMax,\
    minPar,maxEdgeLength,aggOffset,parDiameterList):

    # Store particle nodes that fall inside the bin
    binTestParticles = np.all([(nodes[:,0] > binMin[0]) , \
        (nodes[:,0] < binMax[0]) , (nodes[:,1] > binMin[1]) , \
        (nodes[:,1] < binMax[1]) , (nodes[:,2] > binMin[2]) , \
        (nodes[:,2] < binMax[2])],axis=0)
    existingnodes = nodes[binTestParticles,:]
    existingAggD = parDiameterList[binTestParticles]

    # Compute distance between particles 
    if len(existingnodes>0):
        nodalDistance = np.linalg.norm(center-existingnodes, axis=1)
        aggOffsetDist = nodalDistance - aggDiameter/2 - existingAggD\
            /2 - aggOffset
    else: 
        aggOffsetDist = np.array(([1]))

    # Kill and return if overlap
    if (aggOffsetDist<0).any():
        return True,"NA"

    # Store edge nodes that fall inside the bin
    binTestSurf = np.all([(facePoints[:,0] > binMin[0]) , \
        (facePoints[:,0] < binMax[0]) , (facePoints[:,1] > \
            binMin[1]) , (facePoints[:,1] < binMax[1]) ,\
        (facePoints[:,2] > binMin[2]) , (facePoints[:,2] < \
            binMax[2])],axis=0)     
    existingSurf = facePoints[binTestSurf,:]

    # Compute distance between particle and edge nodes
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
        if (surfNodalDistance<math.sqrt((maxEdgeLength/2)**2+(aggDiameter/2)**2)).any():
            return False,True

    return False,False