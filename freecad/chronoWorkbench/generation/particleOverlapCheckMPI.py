import numpy as np



def overlapCheckMPI(center,aggDiameter,binMin,binMax,minPar,aggOffset,nodes,diameters):

    # Store particle nodes that fall inside the bin
    binTestParticles = np.all([(nodes[:,0] > binMin[0]) , \
        (nodes[:,0] < binMax[0]) , (nodes[:,1] > binMin[1]) , \
        (nodes[:,1] < binMax[1]) , (nodes[:,2] > binMin[2]) , \
        (nodes[:,2] < binMax[2])],axis=0)
    existingnodes = nodes[binTestParticles,:]
    existingAggD = diameters[binTestParticles]

    # Compute distance between particles 
    if len(existingnodes>0):
        nodalDistance = np.linalg.norm(center-existingnodes, axis=1)
        aggOffsetDist = nodalDistance - aggDiameter/2 - existingAggD\
            /2 - aggOffset
    else: 
        aggOffsetDist = np.array(([1]))

    # Kill and return if overlap
    if (aggOffsetDist<0).any():
        return True

    # Otherwise return false
    return False