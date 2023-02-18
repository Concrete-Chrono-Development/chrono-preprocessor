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

import numpy as np

def insideCheck(vertices, center, parDiameter, max_dist):
    distances_squared = np.sum((vertices - center)**2, axis=1)
    if np.all(distances_squared > (parDiameter/2 + max_dist)**2):
        return False
    else:
        return True