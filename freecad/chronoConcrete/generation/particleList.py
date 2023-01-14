import numpy as np



def particleList(parVolTotal,minPar,maxPar,newSieveCurveD,cdf,kappa_i,NewSet,fullerCoef):

    q = 3.0-fullerCoef

    smallParVolume = 4/3*math.pi*(minPar/2)**3
    maxParNum = np.ceil(parVolTotal/smallParVolume)
    parDiameter = np.zeros(int(maxParNum))
    parVol = np.zeros(int(maxParNum))

    i = 0

    # Fuller Coefficient Case
    if newSieveCurveD == 0:

        # Count particles individually
        if len(parDiameter) <= 100:
            while sum(parVol) < parVolTotal:
                parDiameter[i] = minPar*(1-np.random.rand(1)*(1-minPar**q\
                    /maxPar**q))**(-1/q)
                parVol[i] = 4/3*math.pi*(parDiameter[i]/2)**3
                i = i+1         

        # Count particles in bunches of 100
        elif len(parDiameter) <= 1000:
            while sum(parVol) < parVolTotal:
                parDiameter[i:i+100] = minPar*(1-np.random.rand(100)*\
                    (1-minPar**q/maxPar**q))**(-1/q)
                parVol[i:i+100] = 4/3*math.pi*(parDiameter[i:i+100]/2)**3
                i = i+100

        # Count particles in bunches of 1000
        else:
            while sum(parVol) < parVolTotal:
                parDiameter[i:i+1000] = minPar*(1-np.random.rand(1000)*\
                    (1-minPar**q/maxPar**q))**(-1/q)
                parVol[i:i+1000] = 4/3*math.pi*(parDiameter[i:i+1000]/2)**3
                i = i+1000

    # Sieve Curve Case
    else:

        while sum(parVol) < parVolTotal:

            F = np.random.rand(1)

            for x in range(0,NewSet):

                if (F >= cdf[x] and F < cdf[x+1]) :
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

    return maxParNum,parDiameterList