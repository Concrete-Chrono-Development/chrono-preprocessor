


def readInputs(form):

    # Generation Type
    elementType         = "LDPM"

    # Basic Settings
    if form[0].inelasticQuasi.isChecked():
        constitutiveEQ  = "quasiBrittle"
    else:
        constitutiveEQ  = "elastic"
    paramLocation       = form[0].paramLocation.text()
    numCPU              = form[0].numCPUbox.value()
    numIncrements       = form[0].numPIncBox.value()
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

    # Particle Settings
    minPar              = form[2].minPar.value()
    maxPar              = form[2].maxPar.value()        
    fullerCoef          = form[2].fullerCoef.value()  
    sieveCurveDiameter  = form[2].sieveDiameters.text()        
    sieveCurvePassing   = form[2].sievePassing.text()   

    # Mix Design
    wcRatio             = form[3].wcRatio.value()
    densityWater        = form[3].waterDensity.text()
    cementC             = form[3].cementContent.text()
    flyashC             = form[3].flyashContent.text()
    silicaC             = form[3].silicaContent.text()
    scmC                = form[3].scmContent.text()
    cementDensity       = form[3].cementDensity.text()
    flyashDensity       = form[3].flyashDensity.text()
    silicaDensity       = form[3].silicaDensity.text()
    scmDensity          = form[3].scmDensity.text()
    airFrac1            = form[3].airFrac.value()
    airFrac2            = form[3].airFracArb.value()

    # Additional Parameters
    # ... Coming Soon ...


    return elementType, \
        constitutiveEQ, paramLocation, numCPU, numIncrements,placementAlg,\
        geoType, dimensions,\
        minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
        wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
        cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, airFrac2