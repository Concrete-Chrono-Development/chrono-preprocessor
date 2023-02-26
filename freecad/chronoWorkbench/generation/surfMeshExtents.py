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

import numpy as np

def surfMeshExtents(vertices):

    minC = np.amin(vertices, axis=0)
    maxC = np.amax(vertices, axis=0)

    return minC, maxC