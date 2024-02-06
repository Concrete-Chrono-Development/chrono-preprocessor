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

def check_multiMat_matVol(subtetVol,facetMaterial,aggVoxels,itzVoxels,binderVoxels):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - subtetVol:       volume of subtet
    - facetMaterial:   list of all facet materials
    - aggVoxels:       list of all aggregate voxels
    - itzVoxels:       list of all itz voxels
    - binderVoxels:    list of all binder voxels
    --------------------------------------------------------------------------
    ### Outputs ###
    - itzVolFracSim:   simulated ITZ volume fraction
    - binderVolFracSim:simulated binder volume fraction
    - aggVolFracSim:   simulated aggregate volume fraction
    - itzVolFracAct:   actual ITZ volume fraction
    - binderVolFracAct:actual binder volume fraction
    - aggVolFracAct:   actual aggregate volume fraction
    --------------------------------------------------------------------------
    """  

    # Facet volumes
    itzVol = sum(subtetVol*(facetMaterial==1))

    binderVol = sum(subtetVol*(facetMaterial==2))

    aggVol = sum(subtetVol*(facetMaterial==3))

    # Facet volume fractions
    itzVolFracSim = itzVol/(sum((itzVol,binderVol,aggVol)))

    binderVolFracSim = binderVol/(sum((itzVol,binderVol,aggVol)))

    aggVolFracSim = aggVol/(sum((itzVol,binderVol,aggVol)))

    # Voxel volume fractions
    itzVolFracAct = len(itzVoxels)/(len(aggVoxels)+len(itzVoxels)+len(binderVoxels))

    binderVolFracAct = len(binderVoxels)/(len(aggVoxels)+len(itzVoxels)+len(binderVoxels))

    aggVolFracAct = len(aggVoxels)/(len(aggVoxels)+len(itzVoxels)+len(binderVoxels))

    return itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct