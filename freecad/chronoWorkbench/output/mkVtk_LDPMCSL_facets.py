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
## For U.S. Army ERDC Contract No. W9132T22C0015
## Primary Authors: Matthew Troemner
## ===========================================================================
##
## Function to generate VTK file for visualization in Paraview of facets
## and their corresponding material for LDPM models.
##
## ===========================================================================

import numpy as np
from pathlib import Path


def mkVtk_LDPMCSL_facets(geoName,tempPath,facetPointData,facetCellData):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - geoName:          Name of the geometry file
    - tempPath:         Path to the temporary directory
    - facetPointData:   List of facet point data
    - facetCellData:    List of facet cell data
    --------------------------------------------------------------------------
    ### Outputs ###
    - A VTK file that can be visualized in Paraview
    --------------------------------------------------------------------------
    """

    while facetPointData.shape[0] % 3 != 0:
        facetPointData = np.concatenate((facetPointData, np.zeros((1, facetPointData.shape[1]))), axis=0)

    facetPointData = np.around(facetPointData.reshape(-1,9), decimals=6) # reshape and condense to save memory/space

    with open(Path(tempPath + geoName + \
        '-para-facets.000.vtk'),"w") as f:                                                                          
        f.write('# vtk DataFile Version 2.0\n')
        f.write('Facet Visual File\n') 
        f.write('ASCII\n')    
        f.write('DATASET POLYDATA\n')        

        f.write('POINTS ' + str(len(facetPointData)*3) + ' float \n') 
        f.write("\n".join(" ".join(map(str, x)) for x in facetPointData))
        f.write('\n\n')  

        f.write('POLYGONS ' + str(len(facetCellData)) + ' ' \
            + str(round(len(facetCellData)*4)) +'\n3 ')
        f.write("\n3 ".join(" ".join(map(str, x)) for x in facetCellData))

        f.write('\n\n')  

        ## NEED TO FIX BELOW FOR MULTI-MATERIAL IMPLEMENTATION
        #if multiMaterial in ['on','On','Y','y','Yes','yes']:  
        #    f.write('\nCELL_DATA ' + str(len(facetMaterial)) + '\n')
        #    f.write('FIELD FieldData 1\n')
        #    f.write('material 1 ' + str(len(facetMaterial)) + ' float\n')
        #    for x in facetMaterial:
        #        f.write("%s\n" % x)
        #if cementStructure in ['on','On','Y','y','Yes','yes']:  
        #    f.write('\nCELL_DATA ' + str(len(facetMaterial)) + '\n')
        #    f.write('FIELD FieldData 1\n')
        #    f.write('material 1 ' + str(len(facetMaterial)) + ' float\n')
        #    for x in facetMaterial:
        #        f.write("%s\n" % x) 