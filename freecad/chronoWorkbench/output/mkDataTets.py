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
## Function to generate and write a data file of all tetrahedra in an LDPM
## model, for later use in Project Chrono.
##
## ===========================================================================

from pathlib import Path
import numpy as np

def mkDataTets(geoName,tempPath,allTets):
    
    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - geoName:      Name of the geometry file
    - tempPath:     Path to the temporary directory
    - allTets:      Array of all tetrahedra in the model
    --------------------------------------------------------------------------
    ### Outputs ###
    - A data file of all tetrahedra in the model
    --------------------------------------------------------------------------
    """

    np.savetxt(Path(tempPath + geoName + \
        '-data-tets.dat'), allTets.astype(int), fmt='%i', delimiter=' ', comments=''\
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
// Tetrahedra Data File\n\
// ================================================================================\n\
//\n\
// Data Structure:\n\
// Node 1 Node 2 Node 3 Node 4\n\
//\n\
// ================================================================================')


#    with open(Path(tempPath + geoName + \
#        '-data-tets.dat'),"w") as f:                                            
#        f.write('// Mesh Data Generated with LDPM Mesh Generation Tool\n')
#        f.write('\n') 
#        f.write('// Number of Tets: ' + str(len(allTets)) + '\n')  
#        f.write("\n".join(" ".join(map(str, x)) for x in \
#            allTets.astype(int)))   
#        f.write('\n')