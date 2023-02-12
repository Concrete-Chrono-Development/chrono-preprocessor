import numpy as np

def insideCheck(vertices, center, parDiameter, max_dist):
    distances_squared = np.sum((vertices - center)**2, axis=1)
    if np.all(distances_squared > (parDiameter/2 + max_dist)**2):
        return False
    else:
        return True