import numpy as np



def overlapCheckMPI(center,aggDiameter,binMin,binMax,minPar,aggOffset,vertices,diameters):

    # Store particle vertices that fall inside the bin
    binTestParticles = np.all([(vertices[:,0] > binMin[0]) , \
        (vertices[:,0] < binMax[0]) , (vertices[:,1] > binMin[1]) , \
        (vertices[:,1] < binMax[1]) , (vertices[:,2] > binMin[2]) , \
        (vertices[:,2] < binMax[2])],axis=0)
    existingvertices = vertices[binTestParticles,:]
    existingAggD = diameters[binTestParticles]

    # Compute distance between particles 
    if len(existingvertices>0):
        nodalDistance = np.linalg.norm(center-existingvertices, axis=1)
        aggOffsetDist = nodalDistance - aggDiameter/2 - existingAggD\
            /2 - aggOffset
    else: 
        aggOffsetDist = np.array(([1]))

    # Kill and return if overlap
    if (aggOffsetDist<0).any():
        return True

    # Otherwise return false
    return False