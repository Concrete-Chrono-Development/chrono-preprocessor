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


def mkVtk_LDPMCSL_facets(geoName,tempPath,tetFacets):

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


    FacetPoints = tetFacets.reshape(-1,3)

    
    # Make an array like above but goes from 0 to 3xlen(FacetPoints)   
    # This is the cell data for the facets
    FacetCells = np.tile(np.arange(0,len(FacetPoints)),3).reshape(-1,3)





    with open(Path(tempPath + geoName + \
        '-para-facets.000.vtk'),"w") as f:                                                                          
        f.write('# vtk DataFile Version 2.0\n')
        f.write('Facet Visual File\n') 
        f.write('ASCII\n')    
        f.write('DATASET POLYDATA\n')        

        f.write('POINTS ' + str(len(FacetPoints)) + ' float \n') 
        f.write("\n".join(" ".join(map(str, x)) for x in FacetPoints))
        f.write('\n\n')  

        f.write('POLYGONS ' + str(len(FacetCells)) + ' ' \
            + str(round(len(FacetCells)*4)) +'\n3 ')
        f.write("\n3 ".join(" ".join(map(str, x)) for x in FacetCells))

        f.write('\n\n')  