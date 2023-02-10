import numpy as np


def particleVol(wcRatio, volFracAir, fullerCoef, cementC, densityCement, densityWater,
                flyashC, silicaC, scmC, flyashDensity, silicaDensity, scmDensity,
                vertices, tets, minPar, maxPar, sieveCurveDiameter, sieveCurvePassing):
    # Convert density to Kg/m3
    cementC *= 1.0E+12
    densityCement *= 1.0E+12
    densityWater *= 1.0E+12

    # Gets volume of geometry
    tetVolume = meshVolume(vertices, tets)

    # Shift sieve curve if needed
    if sieveCurveDiameter:
        # Shifts sieve curve to appropriate range
        newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max = sieveCurve(minPar, maxPar,
                                                                           sieveCurveDiameter, sieveCurvePassing)
    else:
        newSieveCurveD, newSieveCurveP, w_min, w_max, NewSet = 0, 0, 0, 0, 0

    # Calculates volume of particles needed
    parVolTotal, cdf, cdf1, kappa_i = aggVolume(tetVolume, wcRatio, cementC,
                                                volFracAir, fullerCoef, flyashC, silicaC, scmC, flyashDensity,
                                                silicaDensity, scmDensity,
                                                densityCement, densityWater, minPar, maxPar,
                                                newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max)

    return parVolTotal, cdf, cdf1, kappa_i


def meshVolume(vertices, tets):
    tets = tets.astype(int)

    coord1 = vertices[tets[:, 0] - 1]
    coord2 = vertices[tets[:, 1] - 1]
    coord3 = vertices[tets[:, 2] - 1]
    coord4 = vertices[tets[:, 3] - 1]

    tetVolume = abs(np.vdot(np.transpose(coord1 - coord4),
                            np.transpose(np.cross((coord2 - coord4), (coord3 - coord4)))))/6

    return tetVolume




def sieveCurve(minPar, maxPar, sieveCurveDiameter, sieveCurvePassing):
    nblines = len(sieveCurveDiameter)
    minRange = next(i for i in range(nblines) if minPar >= sieveCurveDiameter[i] and minPar < sieveCurveDiameter[i + 1])
    maxRange = next(i for i in range(nblines) if maxPar > sieveCurveDiameter[i] and maxPar <= sieveCurveDiameter[i + 1])

    w_min = sieveCurvePassing[minRange] + (sieveCurvePassing[minRange + 1] - sieveCurvePassing[minRange]) / \
            (sieveCurveDiameter[minRange + 1] - sieveCurveDiameter[minRange]) * (minPar - sieveCurveDiameter[minRange])
    w_max = sieveCurvePassing[maxRange] + (sieveCurvePassing[maxRange + 1] - sieveCurvePassing[maxRange]) / \
            (sieveCurveDiameter[maxRange + 1] - sieveCurveDiameter[maxRange]) * (maxPar - sieveCurveDiameter[maxRange])
    
    new_set = maxRange - minRange + 1
    newSieveCurveD = [minPar if i == 0 else sieveCurveDiameter[minRange + i] if i > 0 and i < new_set else maxPar for i in range(new_set + 1)]
    newSieveCurveP = [0 if i == 0 else (sieveCurvePassing[minRange + i] - w_min) / (w_max - w_min) if i > 0 and i < new_set else 1 for i in range(new_set + 1)]

    newSieveCurveD = np.trim_zeros(newSieveCurveD, trim='b')
    newSieveCurveP = np.trim_zeros(newSieveCurveP, trim='b')

    return newSieveCurveD, newSieveCurveP, new_set, w_min, w_max




def sieveCurve(minPar, maxPar, sieveCurveDiameter, sieveCurvePassing):
    nblines = len(sieveCurveDiameter)
    
    # Find the indices of the minimum and maximum values in the sieve curve
    minRange = next(i for i in range(nblines) if minPar >= sieveCurveDiameter[i] and minPar < sieveCurveDiameter[i+1])
    maxRange = next(i for i in range(nblines) if maxPar > sieveCurveDiameter[i] and maxPar <= sieveCurveDiameter[i+1])
    
    # Interpolate the minimum and maximum values
    w_min = sieveCurvePassing[minRange] + (sieveCurvePassing[minRange + 1] - sieveCurvePassing[minRange]) / \
            (sieveCurveDiameter[minRange + 1] - sieveCurveDiameter[minRange]) * (minPar - sieveCurveDiameter[minRange])
    w_max = sieveCurvePassing[maxRange] + (sieveCurvePassing[maxRange + 1] - sieveCurvePassing[maxRange]) / \
            (sieveCurveDiameter[maxRange + 1] - sieveCurveDiameter[maxRange]) * (maxPar - sieveCurveDiameter[maxRange])
    
    # Get the number of values in the new sieve curve
    NewSet = maxRange - minRange + 1
    
    # Create arrays for the new sieve curve diameters and passing values
    newSieveCurveD = np.zeros(NewSet + 1)
    newSieveCurveP = np.zeros(NewSet + 1)
    
    # Populate the new sieve curve arrays
    newSieveCurveD[0] = minPar
    newSieveCurveP[0] = 0
    for i in range(1, NewSet):
        newSieveCurveD[i] = sieveCurveDiameter[minRange + i]
        newSieveCurveP[i] = (sieveCurvePassing[minRange + i] - w_min) / (w_max - w_min)
    newSieveCurveD[NewSet] = maxPar
    newSieveCurveP[NewSet] = 1
    
    return newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max



def aggVolume(tetVolume, wcRatio, cementC, volFracAir, fullerCoef, flyashC, silicaC, scmC, flyashDensity, silicaDensity, scmDensity, densityCement, densityWater, minPar, maxPar, newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max):
    # Calculate volume fractions
    volFracCement = cementC / densityCement
    volFracFly = flyashC / flyashDensity
    volFracSilica = silicaC / silicaDensity
    volFracSCM = scmC / scmDensity
    volFracWater = wcRatio * cementC / densityWater
    volFracPar = (1 - volFracCement - volFracFly - volFracSilica - volFracSCM - volFracWater - volFracAir)

    # If newSieveCurveD is 0, calculate volFracParSim using minPar and maxPar
    if newSieveCurveD == 0:
        volFracParSim = (1 - (minPar / maxPar)**fullerCoef) * volFracPar
        kappa_i = []
    else:
        # Calculate Kappa
        A = sum((newSieveCurveP[i+1] - newSieveCurveP[i])/(newSieveCurveD[i+1] - newSieveCurveD[i]) * (1/newSieveCurveD[i]**2 - 1/newSieveCurveD[i+1]**2) for i in range(NewSet))
        kappa = 2 / A
        kappa_i = [kappa * (newSieveCurveP[i+1] - newSieveCurveP[i]) / (newSieveCurveD[i+1] - newSieveCurveD[i]) for i in range(NewSet)]
        kappa_i = np.trim_zeros(kappa_i, trim='b')

        # Calculate volume fraction and cumulative distribution function
        w_sim = w_max - w_min
        volFracParSim = w_sim * volFracPar
        cdf = np.cumsum([kappa_i[i] * (1/newSieveCurveD[i]**2 - 1/newSieveCurveD[i+1]**2)/2 for i in range(NewSet)])
        cdf = np.trim_zeros(cdf, trim='b')

    parVolTotal = volFracParSim * tetVolume

    return parVolTotal, cdf, kappa_i
