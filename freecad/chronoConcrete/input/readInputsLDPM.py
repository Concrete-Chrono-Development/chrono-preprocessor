


def readInputs(form):

    # Generation Type
    elementType         = "LDPM"

    # Basic Settings
    setupFile           = form[0].setupFile.text()
    constitutiveEQ      = form[0].constEQ.currentText()
    matParaSet          = form[0].matParaSet.currentText()

    numCPU              = form[0].numCPUbox.value()
    numIncrements       = form[0].numPIncBox.value()
    maxIter             = form[0].numIncBox.value()
    placementAlg        = form[0].placementAlg.currentText()

    # Geometry Settings
    geoType             = form[1].geometryType.currentText()
    dimensions = []
    if geoType == "Box":
        dimensions.append(form[1].boxLength.text())
        dimensions.append(form[1].boxWidth.text())
        dimensions.append(form[1].boxHeight.text())
    if geoType == "Cylinder":
        dimensions.append(form[1].cylinderHeight.text())
        dimensions.append(form[1].cylinderRadius.text())
    if geoType == "Cone":
        dimensions.append(form[1].coneHeight.text())
        dimensions.append(form[1].coneRadius1.text())
        dimensions.append(form[1].coneRadius2.text())
    if geoType == "Sphere":
        dimensions.append(form[1].sphereRadius.text())
    if geoType == "Ellipsoid":
        dimensions.append(form[1].ellipsoidRadius1.text())
        dimensions.append(form[1].ellipsoidRadius2.text())
        dimensions.append(form[1].ellipsoidRadius3.text())
        dimensions.append(form[1].ellipsoidAngle1.text())
        dimensions.append(form[1].ellipsoidAngle2.text())
        dimensions.append(form[1].ellipsoidAngle3.text())
    if geoType == "Prism":
        dimensions.append(form[1].prismCircumradius.text())
        dimensions.append(form[1].prismHeight.text())
        dimensions.append(form[1].prismPolygon.text())
    if geoType == "Notched Prism":
        dimensions.append(form[1].notchBoxLength.text())
        dimensions.append(form[1].notchBoxWidth.text())
        dimensions.append(form[1].notchBoxHeight.text())
        dimensions.append(form[1].notchWidth.text())
        dimensions.append(form[1].notchDepth.text())
    cadFile             = form[1].cadFile.toPlainText()


    # Particle Settings
    minPar              = float(form[2].minPar.value() or 0)
    maxPar              = float(form[2].maxPar.value() or 0)      
    fullerCoef          = float(form[2].fullerCoef.value() or 0)  
    sieveCurveDiameter  = form[2].sieveDiameters.text()        
    sieveCurvePassing   = form[2].sievePassing.text()   

    # Mix Design
    wcRatio             = float(form[3].wcRatio.value() or 0)
    densityWater        = float(form[3].waterDensity.text() or 0)
    cementC             = float(form[3].cementContent.text() or 0)
    flyashC             = float(form[3].flyashContent.text() or 0)
    silicaC             = float(form[3].silicaContent.text() or 0)
    scmC                = float(form[3].scmContent.text() or 0)
    cementDensity       = float(form[3].cementDensity.text() or 0)
    flyashDensity       = float(form[3].flyashDensity.text() or 0)
    silicaDensity       = float(form[3].silicaDensity.text() or 0)
    scmDensity          = float(form[3].scmDensity.text() or 0)
    airFrac1            = float(form[3].airFrac.value() or 0)
    airFrac2            = float(form[3].airFracArb.value() or 0)

    # Additional Parameters
    # ... Coming Soon ...

    # Generation Data
    outputDir           = form[5].outputDir.text()

    return elementType, \
        setupFile, constitutiveEQ, matParaSet, \
        numCPU, numIncrements,maxIter,placementAlg,\
        geoType, dimensions, cadFile,\
        minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
        wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
        cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, airFrac2,\
        outputDir