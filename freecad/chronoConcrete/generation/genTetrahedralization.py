import os
import numpy as np
from pathlib import Path

from freecad.chronoConcrete import TETGENPATH




def genTetrahedralization(nodes,vertices2D,triangles2D,geoName,verbose,tempPath):
    
    
    

    # Prepare file of internal nodes and external nodes/facets for Tetgen
    nodeRange = np.arange(len(nodes))+1
    nodeList = np.vstack((nodeRange,nodes.T)).T
    
    

    with open(Path(tempPath + geoName + '2D.a.node'),"w") as f:                                       
        f.write(str(len(nodes)) + ' 3 0 0\n ')                                   
        f.write("\n ".join(" ".join(map(str, x)) for x in nodeList))

    print('Starting Tetgen for tetrahedralization construction.')
    
    if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
        tetgenCommand = str(Path(TETGENPATH + '/tetgen')) + ' -pYiO0/1S0k ' \
            + str(Path(tempPath + geoName + '2D.mesh'))
    else:
        tetgenCommand = str(Path(TETGENPATH + '/tetgen')) + ' -pYiO0/1S0kQ ' \
            + str(Path(tempPath + geoName + '2D.mesh'))

    os.system(tetgenCommand)

    try:
        os.rename(Path(tempPath + geoName + '2D.1.vtk'),Path(tempPath + geoName \
            + '-para-mesh.000.vtk'))
    except:
        print("Tetgen failed during tetrahedralization.")
        print("If this issue persists, you may need to use another geometry or particle distribution.")
        




    os.remove(Path(tempPath + geoName + '.mesh'))
    os.remove(Path(tempPath + geoName + '2D.mesh'))
    try:
        os.remove(Path(tempPath + geoName + '2D.stl'))
    except:
        pass
    try:
        os.remove(Path(tempPath + geoName + '2D.1.edge'))
    except:
        pass
    try:
        os.remove(Path(tempPath + geoName + '2D.1.face'))
    except:
        pass
    try:    
        os.remove(Path(tempPath + geoName + '2D.a.node'))
    except:
        pass
    os.rename(Path(tempPath + geoName + '2D.1.ele'),Path(tempPath + geoName \
        + '.ele'))
    os.rename(Path(tempPath + geoName + '2D.1.node'),Path(tempPath + geoName \
        + '.node')) 
    os.replace(Path(tempPath + geoName + '-para-mesh.000.vtk'), \
        Path('meshes/' + geoName + '/' + geoName + '-para-mesh.000.vtk'))
    os.replace(Path(tempPath + geoName + '.ele'), Path('meshes/' + geoName \
        + '/' + geoName + '.ele'))
    os.replace(Path(tempPath + geoName + '.node'), Path('meshes/' + geoName \
        + '/' + geoName + '.node'))