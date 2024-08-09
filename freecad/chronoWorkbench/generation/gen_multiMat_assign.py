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
## Developed by Northwestern University
## For U.S. Army ERDC Contract No. W9132T22C0015
## Primary Authors: Matthew Troemner
## ================================================================================
##
## This file 
##
## ================================================================================

import numpy as np

def gen_multiMat_assign(allNodes, materialList, aggVoxels, itzVoxels, binderVoxels, internalNodes, multiMatX, multiMatY, multiMatZ, multiMatRes, minC):

    # Calculate voxel coordinates for all voxel types at once
    def get_voxel_centers(voxels):
        xVoxels = np.floor(voxels / (multiMatZ * multiMatY))
        yVoxels = np.floor((voxels - xVoxels * multiMatZ * multiMatY) / multiMatZ)
        zVoxels = voxels - xVoxels * multiMatZ * multiMatY - yVoxels * multiMatZ

        voxel_coords = np.stack([xVoxels, yVoxels, zVoxels], axis=-1) * multiMatRes + minC
        voxel_centers = voxel_coords + (multiMatRes / 2) - multiMatRes / 2
        return voxel_centers

    # Get centers for all voxel types
    aggCenters = get_voxel_centers(aggVoxels)
    itzCenters = get_voxel_centers(itzVoxels)
    binderCenters = get_voxel_centers(binderVoxels)

    # Combine centers and materials into one array
    voxel_centers = np.concatenate([aggCenters, itzCenters, binderCenters], axis=0)
    materials = np.concatenate([np.full(len(aggVoxels), 3),
                                np.full(len(itzVoxels), 1),
                                np.full(len(binderVoxels), 2)])

    for i in range(len(allNodes) - len(internalNodes)):
        node = allNodes[i]

        # Calculate distances from the node to all voxel centers at once
        distances = np.linalg.norm(voxel_centers - node, axis=1)

        # Find the index of the closest voxel and assign the corresponding material
        closest_index = np.argmin(distances)
        materialList[i] = materials[closest_index]

    return materialList

