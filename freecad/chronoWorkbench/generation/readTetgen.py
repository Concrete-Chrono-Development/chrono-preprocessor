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


def readTetgen(nodeFile, tetFile):                                       

    allNodes = np.loadtxt(nodeFile, usecols=(1,2,3), \
        skiprows=1)                                   

    allTets = np.loadtxt(tetFile, usecols=(1,2,3,4), \
        skiprows=1)   

    return allNodes, allTets