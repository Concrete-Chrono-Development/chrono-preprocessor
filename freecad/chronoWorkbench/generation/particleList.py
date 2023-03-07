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
## This function generates a sorted array of particle diameters based on user 
## inputs for either a Fuller Coefficient case or a case with a provided
## discrete sieve curve. 
##
## ===========================================================================

import math
import numpy as np


def particleList(parVolTotal, minPar, maxPar, newSieveCurveD, cdf, kappa_i, 
                 NewSet, fullerCoef):
    
    """
    Variable List:
    --------------------------------------------------------------------------
    ### Inputs ###
    parVolTotal:     float, total volume of particles
    minPar:          float, minimum diameter of particles
    maxPar:          float, maximum diameter of particles
    newSieveCurveD:  numpy array, diameters of particles in sieve curve
    cdf:             numpy array, cumulative distribution function
    kappa_i:         numpy array, coefficient used in particle simulation
    NewSet:          int, number of sieves
    fullerCoef:      float, Fuller coefficient
    --------------------------------------------------------------------------
    ### Outputs ###
    maxparNum:       int, maximum number of particles that can be generated
    parDiameterList: numpy array, particle diameters (sorted large-to-small)
    --------------------------------------------------------------------------
    """

    # Determine 'q' value based on Fuller Coefficient
    q = 3.0-fullerCoef

    # Calculate the volume of the smallest particle
    smallparVolume = 4/3*math.pi*(minPar/2)**3
    
    # Determine maximum number of particles
    maxparNum = np.ceil(parVolTotal/smallparVolume)
    
    # Initialize arrays for particle diameters and volumes
    parDiameter = np.zeros(int(maxparNum))
    parVol = np.zeros(int(maxparNum))

    # Counter for particle array indexing
    i = 0

    # If no Sieve Curve is provided
    if newSieveCurveD == 0:

        # Count particles individually for small arrays
        if len(parDiameter) <= 100:
            while sum(parVol) < parVolTotal:
                # Randomly calculate the particle diameter
                parDiameter[i] = minPar*(1-np.random.rand(1)*(1-minPar**q\
                    /maxPar**q))**(-1/q)
                parVol[i] = 4/3*math.pi*(parDiameter[i]/2)**3
                i = i+1         

        # Count particles in bunches of 100 for medium arrays
        elif len(parDiameter) <= 1000:
            while sum(parVol) < parVolTotal:
                # Randomly calculate particle diameters
                parDiameter[i:i+100] = minPar*(1-np.random.rand(100)*\
                    (1-minPar**q/maxPar**q))**(-1/q)
                parVol[i:i+100] = 4/3*math.pi*(parDiameter[i:i+100]/2)**3
                i = i+100

        # Count particles in bunches of 1000 for large arrays
        else:
            while sum(parVol) < parVolTotal:
                # Randomly calculate particle diameters
                parDiameter[i:i+1000] = minPar*(1-np.random.rand(1000)*\
                    (1-minPar**q/maxPar**q))**(-1/q)
                parVol[i:i+1000] = 4/3*math.pi*(parDiameter[i:i+1000]/2)**3
                i = i+1000

    # If a Sieve Curve is provided
    else:
        while sum(parVol) < parVolTotal:
            F = np.random.rand(1)
            for x in range(0,NewSet):
                if (F >= cdf[x] and F < cdf[x+1]) :
                    # Calculate the diameter of the selected particle
                    parDiameter[i] = ((newSieveCurveD[x]**2*kappa_i[x])/(kappa_i[x]-2*(F-cdf[x])*newSieveCurveD[x]**2))**0.5
                    parVol[i] = 4/3*math.pi*(parDiameter[i]/2)**3
                    i = i+1         

    # Remove trailing zeros from arrays
    parDiameter = np.trim_zeros(parDiameter)
    parVol = np.trim_zeros(parVol)

    # Remove accidental extra placed particles
    while sum(parVol) > parVolTotal:
        parDiameter = np.delete(parDiameter, -1)
        parVol = np.delete(parVol, -1)

    # Sort particle diameters large-to-small
    parDiameterList = np.sort(parDiameter)[::-1]

    return maxparNum,parDiameterList