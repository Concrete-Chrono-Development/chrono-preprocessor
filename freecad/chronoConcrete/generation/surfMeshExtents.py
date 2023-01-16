import numpy as np

def surfMeshExtents(vertices):

    minC = np.amin(vertices, axis=0)
    maxC = np.amax(vertices, axis=0)

    return minC, maxC