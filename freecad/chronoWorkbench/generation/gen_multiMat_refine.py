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
##
## ================================================================================


def gen_multiMat_refine(sortedData,itzVolFracSim,\
    binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct,i):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - sortedData:      datastructure with all sorted data
    - itzVolFracSim:   simulated ITZ volume fraction
    - binderVolFracSim:simulated binder volume fraction
    - aggVolFracSim:   simulated aggregate volume fraction
    - itzVolFracAct:   actual ITZ volume fraction
    - binderVolFracAct:actual binder volume fraction
    - aggVolFracAct:   actual aggregate volume fraction
    - i:               index of facet
    --------------------------------------------------------------------------
    ### Outputs ###
    - sortedData:      datastructure with all sorted data
    --------------------------------------------------------------------------
    """  

    # Same material case
    if sortedData[i,3] == sortedData[i,4]:
        #print('0')
        pass

    # Aggregate-ITZ Case
    elif (sortedData[i,3] == 1 and sortedData[i,4] == 3) or \
        (sortedData[i,3] == 3 and sortedData[i,4] == 1):

        if abs(itzVolFracSim-itzVolFracAct) > abs(aggVolFracSim-aggVolFracAct):
            #print('1')
            if itzVolFracSim-itzVolFracAct > 0:
                sortedData[i,2] = 3
            else:
                sortedData[i,2] = 1

        else:
            if aggVolFracSim-aggVolFracAct > 0:
                sortedData[i,2] = 1
            else:
                sortedData[i,2] = 3

    # Aggregate-Binder Case
    elif (sortedData[i,3] == 2 and sortedData[i,4] == 3) or \
        (sortedData[i,3] == 3 and sortedData[i,4] == 2):

        if abs(binderVolFracSim-binderVolFracAct) > abs(aggVolFracSim-aggVolFracAct):
            #print('2')
            if binderVolFracSim-binderVolFracAct > 0:
                sortedData[i,2] = 3
            else:
                sortedData[i,2] = 2

        else:
            if aggVolFracSim-aggVolFracAct > 0:
                sortedData[i,2] = 2
            else:
                sortedData[i,2] = 3

    # ITZ-Binder Case
    elif (sortedData[i,3] == 1 and sortedData[i,4] == 2) or \
        (sortedData[i,3] == 2 and sortedData[i,4] == 1):

        if abs(itzVolFracSim-itzVolFracAct) > abs(binderVolFracSim-binderVolFracAct):
            #print('3')
            if itzVolFracSim-itzVolFracAct > 0:
                sortedData[i,2] = 2
            else:
                sortedData[i,2] = 1

        else:
            if binderVolFracSim-binderVolFracAct > 0:
                sortedData[i,2] = 1
            else:
                sortedData[i,2] = 2
        
    return sortedData