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
## Function to generate a VTK file for visualization in Paraview of a single
## tetrahedron for LDPM and CSL models.
##
## ===========================================================================

from pathlib import Path
import numpy as np



def mkDataParticles(allNodes,parDiameterList,geoName,tempPath):

    # Create diameters list (including zero edge particle diameters)
    allDiameters = np.concatenate((np.array([0.0,]*\
        int(len(allNodes)-len(parDiameterList))),parDiameterList))

    # Generate data file for particles
    with open(Path(tempPath + geoName + \
        '-data-particles.dat'),"w") as f:                                       
        f.write('# Particle Data Generated with LDPM Mesh Generation Tool\n')
        f.write('# [n x y z d]\n')
        f.write('\n')            
        f.write('# Number of Nodes: ' + str(len(allDiameters)) + '\n')    
        f.write('# Number of Aggregates: ' + str(len(parDiameterList)) + '\n')    
        for x in range(0,len(allDiameters)):
            f.write(str(x+1) + ' ' + str(allNodes[x,0]) + ' ' + str(allNodes[x,1]) \
                + ' ' + str(allNodes[x,2]) + ' ' + str(allDiameters[x]) + '\n')
        f.write('\n')  
