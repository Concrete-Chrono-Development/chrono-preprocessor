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
## This functions reads in the CT scan file for fibers
##
## ===========================================================================

import numpy as np

def read_ctScan_file(ctFile):

    """
    Variable List:
    --------------------------------------------------------------------------
    ### Inputs ###
    ctFile:               File path of the CT scan file to read
    --------------------------------------------------------------------------
    ### Outputs ###
    CTScanFiberData:      List of fibers with their properties
    --------------------------------------------------------------------------
    """

    fibersList = np.loadtxt(ctFile, usecols=range(7), \
        skiprows=3)
    CTScanFiberData=[]

    for i in range(len(fibersList)):
        p1Fiber = fibersList[i,1:4]
        p2Fiber = fibersList[i,4:7]

        orienFiber = (p2Fiber-p1Fiber)/\
            np.linalg.norm(p1Fiber-p2Fiber)
        fiberLength = np.linalg.norm(p1Fiber-p2Fiber)

        Data = np.array([p1Fiber[0],p1Fiber[1],p1Fiber[2],p2Fiber[0],p2Fiber[1],p2Fiber[2],\
            orienFiber[0],orienFiber[1],orienFiber[2],fiberLength])
        CTScanFiberData.append(Data)

    CTScanFiberData = np.array(CTScanFiberData).reshape(-1,10)

    return CTScanFiberData    