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
## This file contains the function to check if the voxel system is larger 
## than the geometry.
##
## ===========================================================================

import numpy as np

def sort_multiMat_voxels(multiMatVoxels):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - multiMatVoxels:   Voxel system with multiple materials
    --------------------------------------------------------------------------
    ### Outputs ###
    - aggVoxels:        Voxel numbers for aggregate
    - itzVoxels:        Voxel numbers for ITZ
    - binderVoxels:     Voxel numbers for binder
    --------------------------------------------------------------------------
    """  

    # Number all voxels
    voxelNumbering = np.arange(len(multiMatVoxels))+1

    # Keep only aggregate voxels
    aggFull = (multiMatVoxels>3)*voxelNumbering
    aggVoxels = aggFull[aggFull != 0]
    
    # Keep only ITZ voxels
    itzFull = (multiMatVoxels==2)*voxelNumbering
    itzVoxels = itzFull[itzFull != 0]

    # Keep only binder voxels
    binderFull = (multiMatVoxels==0)*voxelNumbering
    binderVoxels = binderFull[binderFull != 0]

    # Get aggregate voxel IDs
    aggVoxelIDs = multiMatVoxels[aggVoxels-1]

    return aggVoxels,itzVoxels,binderVoxels,aggVoxelIDs