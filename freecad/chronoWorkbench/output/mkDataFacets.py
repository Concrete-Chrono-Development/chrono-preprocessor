## ================================================================================
## CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
##
## Copyright (c) 2023 
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be found
## in the LICENSE file at the top level of the distribution and at
## github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
##
## ================================================================================
## Author: Matthew Troemner
## ================================================================================
##
## Description coming soon...
##
##
## ================================================================================

from pathlib import Path
import numpy as np


def mkDataFacets(geoName,tempPath,facetData,facetPointData):
    

    np.savetxt(Path(tempPath + geoName + \
        '-data-facets.dat'), facetData, fmt='%.10g', comments = '', delimiter=' '\
        ,header='\
// ================================================================================\n\
// CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor\n\
//\n\
// Copyright (c) 2023 \n\
// All rights reserved. \n\
//\n\
// Use of the code that generated this file is governed by a BSD-style license that\n\
// can be found in the LICENSE file at the top level of the distribution and at\n\
// github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE\n\
//\n\
// ================================================================================\n\
// Facet Data File\n\
// ================================================================================\n\
//\n\
// Data Structure:\n\
// Tet IDx IDy IDz Vol pArea cx cy cz px py pz qx qy qz sx sy sz mF\n\
// One line per facet, ordering is Tet 1 (Facet 1-12),...,Tet N (Facet 1-12)\n\
//\n\
// ================================================================================')