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
## This section of code defines several functions for calculating the volume 
## of particles in a mesh, including the meshVolume function for calculating 
## the volume of a mesh composed of tetrahedrons, the sieveCurve function for 
## shifting a total sieve curve to a specific sieve curve, and the parVolume 
## function for calculating the total volume of particles in a mixture and 
## related values. The particleVol function utilizes these functions to 
## calculate the volume of particles in a mesh based on various inputs, such 
## as water to cement ratio, mass of various components, density of those 
## components, and particle size data. 
##
## ===========================================================================

import numpy as np


def particleVol(wcRatio, volFracAir, fullerCoef, cementC, densityCement, 
                densityWater, flyashC, silicaC, scmC, flyashDensity, 
                silicaDensity, scmDensity, fillerC, fillerDensity,
                vertices, tets, minPar, maxPar, 
                sieveCurveDiameter, sieveCurvePassing):
    
    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    wcRatio:            water to cement ratio
    volFracAir:         volume fraction of air in the mix
    fullerCoef:         Fuller coefficient for particle size distribution
    cementC:            cement content
    densityCement:      density of cement
    densityWater:       density of water
    flyashC:            mass of fly ash
    silicaC:            mass of silica fume
    scmC:               mass of supplementary cementitious materials
    fillerC:            mass of arbitrary filler
    flyashDensity:      density of fly ash
    silicaDensity:      density of silica fume
    scmDensity:         density of supplementary cementitious materials
    fillerDensity:      density of arbitrary filler
    vertices:           an array of vertex coordinates for the mesh
    tets:               an array of tetrahedron vertices for the mesh
    minPar:             minimum particle size
    maxPar:             maximum particle size
    sieveCurveDiameter: an array of particle size data
    sieveCurvePassing:  an array of cumulative percent passing data
    --------------------------------------------------------------------------
    ### Outputs ###
    parVolTotal:        volume of particles
    cdf:                cumulative distribution function
    cdf1:               part of the cumulative distribution function
    kappa_i:            kappa coefficient for particle size distribution
    --------------------------------------------------------------------------
    """

    # Convert density to Kg/m3
    cementC = cementC * (1.0E+12)
    flyashC = flyashC * (1.0E+12)
    silicaC = silicaC * (1.0E+12)
    scmC = scmC * (1.0E+12)
    fillerC = fillerC * (1.0E+12)
    densityCement = densityCement * (1.0E+12)
    flyashDensity = flyashDensity * (1.0E+12)
    silicaDensity = silicaDensity * (1.0E+12)
    scmDensity = scmDensity * (1.0E+12)
    fillerDensity = fillerDensity * (1.0E+12)
    densityWater = densityWater * (1.0E+12)

    # Gets volume of geometry
    tetVolume = meshVolume(vertices, tets)

    # Shift sieve curve if needed
    if sieveCurveDiameter != (0 or None or [] or ""):
        # Shifts sieve curve to appropriate range
        [newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max] = sieveCurve(
            minPar, maxPar, sieveCurveDiameter, sieveCurvePassing)
    else:
        newSieveCurveD, newSieveCurveP, w_min, w_max, NewSet = 0, 0, 0, 0, 0

    # Calculates volume of particles needed
    [parVolTotal, cdf, cdf1, kappa_i] = parVolume(tetVolume, wcRatio, cementC,
                                                  volFracAir, fullerCoef, 
                                                  flyashC, silicaC, scmC, fillerC,
                                                  flyashDensity, silicaDensity, 
                                                  scmDensity, fillerDensity, densityCement,
                                                  densityWater, minPar, maxPar,
                                                  newSieveCurveD, newSieveCurveP, 
                                                  NewSet, w_min, w_max)

    return parVolTotal, cdf, cdf1, kappa_i




# Function to calculate the volume of a mesh composed of tetrahedrons
def meshVolume(vertices, tets):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    vertices:          an array of vertex coordinates for the mesh
    tets:              an array of tetrahedron vertices for the mesh
    --------------------------------------------------------------------------
    ### Outputs ###
    sum(tetVolume):    volume of the mesh
    --------------------------------------------------------------------------
    """

    # Convert tets to integer array
    tets = tets.astype(int)
    
    # Get coordinates of vertices for each tetrahedron
    coord1 = vertices[tets[:, 0] - 1]
    coord2 = vertices[tets[:, 1] - 1]
    coord3 = vertices[tets[:, 2] - 1]
    coord4 = vertices[tets[:, 3] - 1]

    # Calculate the volume of all tetrahedrons in one go
    tetVolume = np.abs(np.sum(np.cross((coord2 - coord4),(coord3 - coord4))*(
        coord1 - coord4), axis=1)) / 6

    # Return the sum of all tetrahedron volumes
    return np.sum(tetVolume)





# Function to shift total sieve curve to coarse sieve curve
def sieveCurve(minPar, maxPar, sieveCurveDiameter, sieveCurvePassing):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    minPar:                minimum particle size
    maxPar:                maximum particle size
    sieveCurveDiameter:    an array of particle size data
    sieveCurvePassing:     an array of cumulative percent passing data
    --------------------------------------------------------------------------
    ### Outputs ###
    newSieveCurveD:        shifted sieve curve diameter data
    newSieveCurveP:        shifted sieve curve percent passing data
    NewSet:                number of lines in shifted sieve curve data
    w_min:                 minimum weight of particles
    w_max:                 maximum weight of particles
    --------------------------------------------------------------------------
    """

    
    # Get number of lines in sieve curve data
    nblines = len(sieveCurveDiameter)
    
    # Find range of diameters containing the minimum particle size
    for i in range(0, nblines):
        if minPar >= sieveCurveDiameter[i] and minPar<sieveCurveDiameter[i+1]:
            minRange = i
            
    # Find range of diameters containing the maximum particle size
    for i in range(0, nblines):
        if maxPar > sieveCurveDiameter[i] and maxPar<=sieveCurveDiameter[i+1]:
            maxRange = i
            
    # Calculate minimum and maximum weights from sieve curve data
    w_min = sieveCurvePassing[minRange] + (
        sieveCurvePassing[minRange+1] - sieveCurvePassing[minRange]) / (
        sieveCurveDiameter[minRange+1] - sieveCurveDiameter[minRange]) * (
        minPar - sieveCurveDiameter[minRange])
    w_max = sieveCurvePassing[maxRange] + (
        sieveCurvePassing[maxRange+1] - sieveCurvePassing[maxRange]) / (
        sieveCurveDiameter[maxRange+1] - sieveCurveDiameter[maxRange]) * (
        maxPar - sieveCurveDiameter[maxRange])
    
    # Get number of data points in new sieve curve
    NewSet = maxRange - minRange + 1
    
    # Initialize new sieve curve arrays
    newSieveCurveD = [0] * 100
    newSieveCurveP = [0] * 100
    
    # Fill in new sieve curve data
    for i in range(0, NewSet+1):
        if i == 0:
            newSieveCurveD[i] = minPar
            newSieveCurveP[i] = 0
        elif NewSet > 1 and i > 0 and i < NewSet:
            newSieveCurveD[i] = sieveCurveDiameter[minRange + i]
            newSieveCurveP[i] = (sieveCurvePassing[minRange + i] - w_min) / (
                w_max - w_min)
        elif NewSet == i:
            newSieveCurveD[i] = maxPar
            newSieveCurveP[i] = 1
            
    # Trim excess zeros from new sieve curve arrays
    newSieveCurveD = np.trim_zeros(newSieveCurveD, trim='b')
    newSieveCurveP = np.trim_zeros(newSieveCurveP, trim='b')
    
    # Return new sieve curve data and related values
    return newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max




# Function to calculate the total volume of particles in a mixture and related
# values
def parVolume(tetVolume, wcRatio, cementC, volFracAir, fullerCoef, flyashC,
              silicaC, scmC, fillerC, flyashDensity, silicaDensity, scmDensity,
              fillerDensity, densityCement, densityWater, minPar, maxPar, 
              newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max):
    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    tetVolume:      total volume of the mix
    wcRatio:        water to cement ratio
    cementC:        cement content
    volFracAir:     volume fraction of air in the mix
    fullerCoef:     Fuller coefficient for particle size distribution
    flyashC:        mass of fly ash
    silicaC:        mass of silica fume
    scmC:           mass of supplementary cementitious materials
    fillerC:        mass of arbitrary filler
    flyashDensity:  density of fly ash
    silicaDensity:  density of silica fume
    scmDensity:     density of supplementary cementitious materials
    fillerDensity:  density of arbitrary filler
    densityCement:  density of cement
    densityWater:   density of water
    minPar:         minimum particle size
    maxPar:         maximum particle size
    newSieveCurveD: an array containing particle size data
    newSieveCurveP: an array containing cumulative percent passing data
    NewSet:         number of data points in newSieveCurveD and newSieveCurveP
    w_min:          minimum value for a weight distribution
    w_max:          maximum value for a weight distribution
    --------------------------------------------------------------------------
    ### Outputs ###
    parVolSimTotal: total volume of particles in the mix
    cdf:            cumulative distribution function
    cdf1:           cumulative distribution function
    kappa_i:        kappa value
    --------------------------------------------------------------------------
    """
    
    # Calculate volume fractions of each component in the mixture and make sure to ignore if density is zero
    if densityCement != 0:
        volFracCement = cementC / densityCement
    else:
        volFracCement = 0
    if flyashDensity != 0:
        volFracFly = flyashC / flyashDensity
    else:
        volFracFly = 0
    if silicaDensity != 0:
        volFracSilica = silicaC / silicaDensity
    else:
        volFracSilica = 0
    if scmDensity != 0:
        volFracSCM = scmC / scmDensity
    else:
        volFracSCM = 0
    if fillerDensity != 0:
        volFracFiller = fillerC / fillerDensity
    else:
        volFracFiller = 0
    if densityWater != 0:
        volFracWater = wcRatio * cementC / densityWater
    else:
        volFracWater = 0

    # Calculate volume fraction of particles
    volFracPar = (1-volFracCement-volFracFly-volFracSilica-volFracSCM-
                  volFracFiller-volFracWater - volFracAir)
    
    # Check if a pre-defined particle size distribution curve is available
    if newSieveCurveD == 0:
        
        # Calculate volume fraction of particles using Fuller's equation
        volFracParSim = (1 - (minPar / maxPar) ** fullerCoef) * volFracPar
        cdf = 0
        cdf1 = 0
        kappa_i = 0
        
    else:
        
        # Calculate Kappa value
        A = 0
        for i in range(0, NewSet):
            B = (newSieveCurveP[i+1] - newSieveCurveP[i]) / (
                newSieveCurveD[i+1] - newSieveCurveD[i]) * (
                1 / newSieveCurveD[i] / newSieveCurveD[i] -
                1 / newSieveCurveD[i+1] / newSieveCurveD[i+1])
            A = A + B
        kappa = 2 / A
        kappa_i = [0] * 100
        for i in range(0, NewSet):
            kappa_i[i] = (kappa * (newSieveCurveP[i+1] - newSieveCurveP[i]) /
                          (newSieveCurveD[i+1] - newSieveCurveD[i]))
        kappa_i = np.trim_zeros(kappa_i, trim='b')
        
        # Calculate volume fraction of particles using simulated values
        w_sim = w_max - w_min
        volFracParSim = w_sim * volFracPar
        
        # Calculate CDF
        cdf = [0] * 100
        cdf1 = [0] * 100
        for i in range(0, NewSet):
            cdf1[i] = (kappa_i[i] * (1 / (newSieveCurveD[i] ** 2) -
                                      1 / (newSieveCurveD[i+1] ** 2))) / 2
            if i > 0:
                cdf[i] = cdf1[i-1] + cdf[i-1]
        cdf[NewSet] = cdf1[NewSet-1] + cdf[NewSet-1]
        cdf = np.trim_zeros(cdf, trim='b')
        cdf1 = np.trim_zeros(cdf1, trim='b')
        
    # Calculate total volume of particles
    parVolSimTotal = volFracParSim * tetVolume

    # Return total particle volume, CDF, and other related values
    return parVolSimTotal, cdf, cdf1, kappa_i














