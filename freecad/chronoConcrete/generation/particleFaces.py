
import numpy as np




# Pulls the coordinates of each external triangle in the mesh 
def particleFaces(vertices,faces):

    faces = faces.astype(int)
    
    coord1 = vertices[faces[:,0]-1]
    coord2 = vertices[faces[:,1]-1]
    coord3 = vertices[faces[:,2]-1] 

    facePoints = np.concatenate((coord1,coord2,coord3))

    return facePoints