## ===========================================================================
## CHRONO WORKBENCH:github.com/Concrete-Chrono-Development/chrono-preprocessor
##
## Copyright (c) 2023 
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be
## found in the LICENSE file at the top level of the distribution and at
## github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
##
## ===========================================================================
## Developed by Northwestern University
## Primary Authors: Matthew Troemner
## ===========================================================================
##
## Description coming soon...
##
##
## ===========================================================================

import os
import numpy as np
from pathlib import Path

from freecad.chronoWorkbench import TETGENPATH




def genTetrahedralization(nodes,vertices2D,triangles2D,geoName,verbose,tempPath):
    
    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - nodes:            Internal nodes of the geometry
    - vertices2D:       External nodes of the geometry
    - triangles2D:      External faces of the geometry
    - geoName:          Name of the geometry
    - verbose:          Print additional information
    - tempPath:         Path to the temporary folder
    --------------------------------------------------------------------------
    ### Outputs ###
    - None
    --------------------------------------------------------------------------
    """  

    # Prepare file of internal nodes and external nodes/facets for Tetgen
    nodeRange = np.arange(len(nodes))+1
    nodeList = np.vstack((nodeRange,nodes.T)).T
    


    # Make external triangles file for Tetgen
    with open(Path(tempPath + geoName + '2D.mesh'),"w") as f:                                       
        f.write('MeshVersionFormatted 2\n')   
        f.write('Dimension\n')   
        f.write('3\n')   
        f.write('Vertices\n')   
        f.write(str(len(vertices2D)) + '\n')                                   
        for x in range(0,len(vertices2D)):
                f.write(str(vertices2D[x,0]) + '    ' \
                    + str(vertices2D[x,1]) + '    '  + str(vertices2D[x,2]) \
                    + '    0\n')
        f.write('Triangles\n')
        f.write(str(len(triangles2D)) + '\n')  
        for x in range(0,len(triangles2D)):
                f.write(str(triangles2D[x,0]) + '    ' \
                    + str(triangles2D[x,1]) + '    '  + str(triangles2D[x,2]) \
                    + '\n')
        f.write('End\n')

    
    # Make internal node file for Tetgen
    with open(Path(tempPath + geoName + '2D.a.node'),"w") as f:                                       
        f.write(str(len(nodes)) + ' 3 0 0\n ')                                   
        f.write("\n ".join(" ".join(map(str, x)) for x in nodeList))




    # Run Tetgen with appropriate switches
    tetgenCommand = str(Path(TETGENPATH + '/tetgen')) + ' -pYiO0/1S0kQ ' \
        + str(Path(tempPath + geoName + '2D.mesh'))
    os.system(tetgenCommand)

    # Check if Tetgen ran by trying to rename the output VTK
    try:
        os.rename(Path(tempPath + geoName + '2D.1.vtk'),Path(tempPath + geoName \
            + '-para-mesh.vtk'))
    except:
        print("Tetgen failed during tetrahedralization.")
        print("If this issue persists, you may need to use another geometry or particle distribution.")
        
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