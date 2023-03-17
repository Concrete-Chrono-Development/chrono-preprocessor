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
## Function to generate VTK file for visualization in Paraview of facets for a
## single polyhedral cell in LDPM and CSL.
##
## ===========================================================================

import numpy as np
from pathlib import Path


def mkVtkSingleCell(allNodes,allTets,parDiameterList,tetFacets,geoName,tempPath):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - allNodes:     List of all nodes
    - allTets:      List of all tetrahedra
    - parDiameterList:  List of particle diameters
    - tetFacets:    List of facets for each tetrahedron
    - geoName:      Name of the geometry
    - tempPath:     Path to the temporary directory
    --------------------------------------------------------------------------
    ### Outputs ###
    - A VTK file that can be visualized in Paraview
    --------------------------------------------------------------------------
    """

    # Select a random non-surface particle from particle diameter list
    randIndex = np.random.randint(int(len(allNodes)-len(parDiameterList)),len(allNodes))
    print(randIndex)


    # Find locations in tet list which contain the given vertex
    location = np.argwhere(allTets==randIndex)   
    print(location)

    # Initialize empty facet list for vertex
    cellFacetNodes=np.zeros(([len(location)*6,9]))

    # Store facets for given cell. Using facet indexing from genTessellationLDPM.py (same as Cusatis 2011a)
    for x in range(0,len(location)):
        if location[x,1] == 0:
            facets = [0,1,2,3,10,11]
        elif location[x,1] == 1:
            facets = [0,1,4,5,6,7]
        elif location[x,1] == 2:
            facets = [2,3,4,5,8,9]
        elif location[x,1] == 3:
            facets = [6,7,8,9,10,11]
        print(facets)

        # Stores facets. One line per facet (9 values for 3 nodes) 
        cellFacetNodes[x*6:x*6+6,:]     = tetFacets[location[x,0],facets[0]*9:facets[0]*9+9]
        cellFacetNodes[x*6+1:x*6+6+1,:] = tetFacets[location[x,0],facets[1]*9:facets[1]*9+9]
        cellFacetNodes[x*6+2:x*6+6+2,:] = tetFacets[location[x,0],facets[2]*9:facets[2]*9+9]
        cellFacetNodes[x*6+3:x*6+6+3,:] = tetFacets[location[x,0],facets[3]*9:facets[3]*9+9]
        cellFacetNodes[x*6+4:x*6+6+4,:] = tetFacets[location[x,0],facets[4]*9:facets[4]*9+9]
        cellFacetNodes[x*6+5:x*6+6+5,:] = tetFacets[location[x,0],facets[5]*9:facets[5]*9+9]
    print(np.shape(tetFacets))  

    # Make triangle connectivity for facets
    cellFacets = np.arange(len(cellFacetNodes)*3).reshape(-1,3).astype(int)

    with open(Path(tempPath + geoName + \
        '-para-singleCellFacets.000.vtk'),"w") as f:                                                                          
        f.write('# vtk DataFile Version 2.0\n')
        f.write('Facet Visual File\n') 
        f.write('ASCII\n')    
        f.write('DATASET POLYDATA\n')        

        f.write('POINTS ' + str(int(len(cellFacetNodes)*3)) + ' float \n') 
        f.write("\n".join(" ".join(map(str, x)) for x in cellFacetNodes))
        f.write('\n\n')  

        f.write('POLYGONS ' + str(len(cellFacets)) + ' ' \
            + str(round(len(cellFacets)*4)) +'\n3 ')
        f.write("\n3 ".join(" ".join(map(str, x)) for x in cellFacets))

        f.write('\n\n')   
