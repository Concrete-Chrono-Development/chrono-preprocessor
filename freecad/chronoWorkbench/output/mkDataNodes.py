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
## Function to generate and write a data file of all nodes in an LDPM model, 
## for later use in Project Chrono.
##
## ===========================================================================

from pathlib import Path
import numpy as np


def mkDataNodes(geoName,tempPath,allNodes):
    
    np.savetxt(Path(tempPath + geoName + \
        '-data-nodes.dat'), allNodes, fmt='%.10g', delimiter=' ', comments=''\
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
// Node Data File\n\
// ================================================================================\n\
//\n\
// Data Structure:\n\
// X Y Z \n\
//\n\
// ================================================================================')