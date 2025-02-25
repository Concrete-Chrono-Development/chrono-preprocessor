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
##
##
## ===========================================================================

import os
from pathlib import Path

# Importing: input
from freecad.chronoWorkbench.input.read_LDPMCSL_inputs                    import read_LDPMCSL_inputs
from freecad.chronoWorkbench.input.read_SPHDEM_inputs                     import read_SPHDEM_inputs

def mkParameters(self,elementSet,tempPath):

    # Read in inputs from input panel
    if elementSet == "LDPMCSL":
        [setupFile, constitutiveEQ, matParaSet, \
            numCPU, numIncrements,maxIter,placementAlg,\
            geoType, dimensions, cadFile,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
            cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, \
            fillerC, fillerDensity, airFrac2,\
            htcToggle, htcLength,\
            fiberToggle, fiberCutting, fiberDiameter, fiberLength, fiberVol, fiberOrientation1, fiberOrientation2, fiberOrientation3, fiberPref, fiberFile, fiberIntersections,\
            multiMatToggle,aggFile,multiMatFile,multiMatRule,\
            grainAggMin, grainAggMax, grainAggFuller, grainAggSieveD, grainAggSieveP,\
            grainITZMin, grainITZMax, grainITZFuller, grainITZSieveD, grainITZSieveP,\
            grainBinderMin, grainBinderMax, grainBinderFuller, grainBinderSieveD, grainBinderSieveP,\
            periodicToggle,\
            outDir, dataFilesGen, visFilesGen, singleTetGen, modelType] = read_LDPMCSL_inputs(self.form)
    else:
        [setupFile, \
            numCPU, numIncrements,maxIter,placementAlg,\
            geoType, dimensions, cadFile,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing, minDistCoef,\
            wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
            cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, \
            fillerC, fillerDensity, airFrac2,\
            outDir, modelType] = read_SPHDEM_inputs(self.form)


    # Make output directory if does not exist
    try:
        os.mkdir(outDir)
    except:
        pass

    # If tempPath is "writeOnly" (meaning we only write the parameter file and don't generate the model) then write to default location
    if tempPath == "writeOnly":
        usePath = Path(outDir + "/chronoWorkbench.cwPar")
    else:
        usePath = Path(tempPath + "/chronoWorkbench.cwPar")

    # Write parameters to file
    if elementSet == "LDPMCSL":
        with open(usePath, "w") as f:
            f.write("""
// ================================================================================
// CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
//
// Copyright (c) 2023 
// All rights reserved. 
//
// Use of the code that generated this file is governed by a BSD-style license that
// can be found in the LICENSE file at the top level of the distribution and at
// github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
//
// ================================================================================
// Chrono Workbench Parameter File
// ================================================================================
// 
// Chrono Workbench developed by Northwestern University
//
// ================================================================================
        \n\n""")
            f.write("constitutiveEQ = " + constitutiveEQ + "\n")
            f.write("matParaSet = " + matParaSet + "\n")
            f.write("numCPU = " + str(numCPU) + "\n")
            f.write("numIncrements = " + str(numIncrements) + "\n")
            f.write("maxIter = " + str(maxIter) + "\n")
            f.write("placementAlg = " + placementAlg + "\n")
            f.write("geoType = " + geoType + "\n")
            f.write("dimensions = " + str(dimensions) + "\n")
            f.write("cadFile = " + cadFile + "\n")
            f.write("minPar = " + str(minPar) + "\n")
            f.write("maxPar = " + str(maxPar) + "\n")
            f.write("fullerCoef = " + str(fullerCoef) + "\n")
            f.write("sieveCurveDiameter = " + str(sieveCurveDiameter) + "\n")
            f.write("sieveCurvePassing = " + str(sieveCurvePassing) + "\n")
            f.write("wcRatio = " + str(wcRatio) + "\n")
            f.write("densityWater = " + str(densityWater) + "\n")
            f.write("cementC = " + str(cementC) + "\n")
            f.write("flyashC = " + str(flyashC) + "\n")
            f.write("silicaC = " + str(silicaC) + "\n")
            f.write("scmC = " + str(scmC) + "\n")
            f.write("cementDensity = " + str(cementDensity) + "\n")
            f.write("flyashDensity = " + str(flyashDensity) + "\n")
            f.write("silicaDensity = " + str(silicaDensity) + "\n")
            f.write("scmDensity = " + str(scmDensity) + "\n")
            f.write("airFrac1 = " + str(airFrac1) + "\n")
            f.write("fillerC = " + str(fillerC) + "\n")
            f.write("fillerDensity = " + str(fillerDensity) + "\n")
            f.write("airFrac2 = " + str(airFrac2) + "\n")
            f.write("htcToggle = " + htcToggle + "\n")
            f.write("htcLength = " + str(htcLength) + "\n")
            f.write("fiberToggle = " + fiberToggle + "\n")
            f.write("fiberCutting = " + fiberCutting + "\n")
            f.write("fiberDiameter = " + str(fiberDiameter) + "\n")
            f.write("fiberLength = " + str(fiberLength) + "\n")
            f.write("fiberVol = " + str(fiberVol) + "\n")
            f.write("fiberOrientation1 = " + str(fiberOrientation1) + "\n")
            f.write("fiberOrientation2 = " + str(fiberOrientation2) + "\n")
            f.write("fiberOrientation3 = " + str(fiberOrientation3) + "\n")
            f.write("fiberPref = " + str(fiberPref) + "\n")
            f.write("fiberFile = " + fiberFile + "\n")
            f.write("multiMatToggle = " + multiMatToggle + "\n")
            f.write("multiMatFile = " + multiMatFile + "\n")
            f.write("aggFile = " + aggFile + "\n")
            f.write("multiMatRule = " + str(multiMatRule) + "\n")
            f.write("grainAggMin = " + str(grainAggMin) + "\n")
            f.write("grainAggMax = " + str(grainAggMax) + "\n")
            f.write("grainAggFuller = " + str(grainAggFuller) + "\n")
            f.write("grainAggSieveD = " + str(grainAggSieveD) + "\n")
            f.write("grainAggSieveP = " + str(grainAggSieveP) + "\n")
            f.write("grainITZMin = " + str(grainITZMin) + "\n")
            f.write("grainITZMax = " + str(grainITZMax) + "\n")
            f.write("grainITZFuller = " + str(grainITZFuller) + "\n")
            f.write("grainITZSieveD = " + str(grainITZSieveD) + "\n")
            f.write("grainITZSieveP = " + str(grainITZSieveP) + "\n")
            f.write("grainBinderMin = " + str(grainBinderMin) + "\n")
            f.write("grainBinderMax = " + str(grainBinderMax) + "\n")
            f.write("grainBinderFuller = " + str(grainBinderFuller) + "\n")
            f.write("grainBinderSieveD = " + str(grainBinderSieveD) + "\n")
            f.write("grainBinderSieveP = " + str(grainBinderSieveP) + "\n")
            f.write("periodicToggle = " + periodicToggle + "\n")
            f.write("dataFilesGen = " + str(dataFilesGen) + "\n")
            f.write("visFilesGen = " + str(visFilesGen) + "\n")
            f.write("singleTetGen = " + str(singleTetGen) + "\n")
            f.write("modelType = " + modelType + "\n")
            f.write("outputDir = " + outDir + "\n")
        print("Parameters written to file")


    elif elementSet == "SPHDEM":
        with open(usePath, "w") as f:
            f.write("""
// ================================================================================
// CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
//
// Copyright (c) 2023 
// All rights reserved. 
//
// Use of the code that generated this file is governed by a BSD-style license that
// can be found in the LICENSE file at the top level of the distribution and at
// github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
//
// ================================================================================
// Chrono Workbench Parameter File
// ================================================================================
// 
// Chrono Workbench developed by Northwestern University
//
// ================================================================================
        \n\n""")
            f.write("numCPU = " + str(numCPU) + "\n")
            f.write("numIncrements = " + str(numIncrements) + "\n")
            f.write("maxIter = " + str(maxIter) + "\n")
            f.write("placementAlg = " + placementAlg + "\n")
            f.write("geoType = " + geoType + "\n")
            f.write("dimensions = " + str(dimensions) + "\n")
            f.write("cadFile = " + cadFile + "\n")
            f.write("minPar = " + str(minPar) + "\n")
            f.write("maxPar = " + str(maxPar) + "\n")
            f.write("fullerCoef = " + str(fullerCoef) + "\n")
            f.write("sieveCurveDiameter = " + str(sieveCurveDiameter) + "\n")
            f.write("sieveCurvePassing = " + str(sieveCurvePassing) + "\n")
            f.write("minDistCoef = " + str(minDistCoef) + "\n")
            f.write("wcRatio = " + str(wcRatio) + "\n")
            f.write("densityWater = " + str(densityWater) + "\n")
            f.write("cementC = " + str(cementC) + "\n")
            f.write("flyashC = " + str(flyashC) + "\n")
            f.write("silicaC = " + str(silicaC) + "\n")
            f.write("scmC = " + str(scmC) + "\n")
            f.write("cementDensity = " + str(cementDensity) + "\n")
            f.write("flyashDensity = " + str(flyashDensity) + "\n")
            f.write("silicaDensity = " + str(silicaDensity) + "\n")
            f.write("scmDensity = " + str(scmDensity) + "\n")
            f.write("airFrac1 = " + str(airFrac1) + "\n")
            f.write("fillerC = " + str(fillerC) + "\n")
            f.write("fillerDensity = " + str(fillerDensity) + "\n")
            f.write("airFrac2 = " + str(airFrac2) + "\n")
            f.write("modelType = " + modelType + "\n")
            f.write("outputDir = " + outDir + "\n")

        print("Parameters written to file")
    
    else:
        print("Error: Element set not recognized")