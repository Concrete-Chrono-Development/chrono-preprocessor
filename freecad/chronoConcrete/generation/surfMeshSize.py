import numpy as np


def surfMeshSize(vertices,faces):

    # Calculate surface mesh max edge length
    maxEdgeLength1  = max(np.linalg.norm(vertices[faces[:,1].astype(int)-1,0:3]-vertices[faces[:,0].astype(int)-1,0:3], axis=1))
    maxEdgeLength2  = max(np.linalg.norm(vertices[faces[:,2].astype(int)-1,0:3]-vertices[faces[:,1].astype(int)-1,0:3], axis=1))
    maxEdgeLength3  = max(np.linalg.norm(vertices[faces[:,2].astype(int)-1,0:3]-vertices[faces[:,0].astype(int)-1,0:3], axis=1))
    maxEdgeLength   = max([maxEdgeLength1,maxEdgeLength2,maxEdgeLength3])

    return maxEdgeLength
