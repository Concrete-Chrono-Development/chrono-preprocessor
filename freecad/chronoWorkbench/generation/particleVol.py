import numpy as np



def particleVol(wcRatio,volFracAir,fullerCoef,cementC,densityCement,densityWater,\
    flyashC,silicaC,scmC,flyashDensity,silicaDensity,scmDensity,\
	vertices,tets,minPar,maxPar,sieveCurveDiameter,sieveCurvePassing):

    # Convert density to Kg/m3
    cementC = cementC*(1.0E+12)
    densityCement = densityCement*(1.0E+12)
    densityWater = densityWater*(1.0E+12)


    # Gets volume of geometry
    tetVolume = meshVolume(vertices,tets)
    

    # Shift sieve curve if needed
    if sieveCurveDiameter != (0 or None or [] or ""):
        # Shifts sieve curve to appropriate range
        [newSieveCurveD,newSieveCurveP,NewSet,w_min,w_max] = sieveCurve(minPar,maxPar,\
            sieveCurveDiameter,sieveCurvePassing)

    else:
        newSieveCurveD, newSieveCurveP, w_min, w_max, NewSet = 0, 0, 0, 0, 0



    # Calculates volume of particles needed
    [parVolTotal,cdf,cdf1,kappa_i] = aggVolume(tetVolume,wcRatio,cementC,\
        volFracAir,fullerCoef,flyashC,silicaC,scmC,flyashDensity,\
        silicaDensity,scmDensity,\
        densityCement,densityWater,minPar,maxPar,\
        newSieveCurveD,newSieveCurveP,NewSet,w_min,w_max)


    return parVolTotal,cdf,cdf1,kappa_i






def meshVolume(vertices,tets):
    tets = tets.astype(int)
    coord1 = vertices[tets[:,0]-1]
    coord2 = vertices[tets[:,1]-1]
    coord3 = vertices[tets[:,2]-1]
    coord4 = vertices[tets[:,3]-1]

    # Calculate the volume of all tetrahedrons in one go
    tetVolume = np.abs(np.sum(np.cross((coord2-coord4),(coord3-coord4)) * (coord1-coord4), axis=1))/6

    return np.sum(tetVolume)



# Shift total sieve curve to coarse sieve curve
def sieveCurve(minPar,maxPar,sieveCurveDiameter,sieveCurvePassing):
    
    nblines = len(sieveCurveDiameter)

    for i in range(0,nblines):
        if minPar >= sieveCurveDiameter[i] and minPar < sieveCurveDiameter[i+1]:      
            minRange = i

    for i in range(0,nblines):
        if maxPar > sieveCurveDiameter[i] and maxPar <= sieveCurveDiameter[i+1]:              
            maxRange = i

    w_min = sieveCurvePassing[minRange]+(sieveCurvePassing[minRange+1]-sieveCurvePassing[minRange])/\
        (sieveCurveDiameter[minRange+1]-sieveCurveDiameter[minRange])*(minPar-sieveCurveDiameter[minRange])
    w_max = sieveCurvePassing[maxRange]+(sieveCurvePassing[maxRange+1]-sieveCurvePassing[maxRange])/\
        (sieveCurveDiameter[maxRange+1]-sieveCurveDiameter[maxRange])*(maxPar-sieveCurveDiameter[maxRange])

    NewSet = maxRange-minRange+1

    newSieveCurveD = [0]*100
    newSieveCurveP = [0]*100

    for i in range(0,NewSet+1):
        if i == 0:
            newSieveCurveD[i] = minPar
            newSieveCurveP[i] = 0
        elif NewSet > 1 and i > 0 and i < NewSet:
            newSieveCurveD[i] = sieveCurveDiameter[minRange + i]
            newSieveCurveP[i] = (sieveCurvePassing[minRange + i] - w_min)/(w_max - w_min)
        elif NewSet == i:
            newSieveCurveD[i] = maxPar
            newSieveCurveP[i] = 1

    newSieveCurveD = np.trim_zeros(newSieveCurveD,trim='b')
    newSieveCurveP = np.trim_zeros(newSieveCurveP,trim='b')

    return newSieveCurveD,newSieveCurveP,NewSet,w_min,w_max


def aggVolume(tetVolume,wcRatio,cementC,volFracAir,fullerCoef,\
    flyashC,silicaC,scmC,flyashDensity,silicaDensity,scmDensity,\
    densityCement,densityWater,minPar,maxPar,newSieveCurveD,newSieveCurveP,NewSet,w_min,w_max):



    volFracCement = cementC/densityCement
    volFracFly = flyashC/flyashDensity
    volFracSilica = silicaC/silicaDensity
    volFracSCM = scmC/scmDensity
    volFracWater = wcRatio*cementC/densityWater
    volFracPar = (1-volFracCement-volFracFly-volFracSilica-volFracSCM-volFracWater-volFracAir)


    if newSieveCurveD == 0:

        volFracParSim = (1-(minPar/maxPar)**fullerCoef)*volFracPar

        cdf = 0
        cdf1 = 0
        kappa_i = 0

    else:

        # Calculate Kappa
        A = 0 

        for i in range(0,NewSet):
            B = (newSieveCurveP[i+1] - newSieveCurveP[i])/(newSieveCurveD[i+1] - newSieveCurveD[i])\
                *(1/newSieveCurveD[i]/newSieveCurveD[i] - 1/newSieveCurveD[i+1]/newSieveCurveD[i+1])
            A = A+B

        kappa = 2/A
        kappa_i = [0]*100

        for i in range(0,NewSet):
            kappa_i[i] = (kappa * (newSieveCurveP[i+1] - newSieveCurveP[i]) / (newSieveCurveD[i+1] - newSieveCurveD[i]))
        kappa_i = np.trim_zeros(kappa_i,trim='b')

        # Calculate Volume Fraction
        w_sim = w_max-w_min
        volFracParSim = w_sim*volFracPar

        # Get CDF
        cdf = [0]*100
        cdf1 = [0]*100
        for i in range(0,NewSet):
            cdf1[i] = kappa_i[i] * (1./(newSieveCurveD[i]**2) - 1/(newSieveCurveD[i+1]**2))/2
            if i > 0:
                cdf[i] = cdf1[i-1] + cdf[i-1]
        cdf[NewSet] = cdf1[NewSet-1] + cdf[NewSet-1]

        cdf = np.trim_zeros(cdf,trim='b')
        cdf1 = np.trim_zeros(cdf1,trim='b')

    parVolTotal = volFracParSim*tetVolume

    return parVolTotal,cdf,cdf1,kappa_i














