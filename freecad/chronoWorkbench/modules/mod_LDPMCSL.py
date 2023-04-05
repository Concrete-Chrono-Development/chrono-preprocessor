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
## Description coming soon...
##
##
## ===========================================================================

# pyright: reportMissingImports=false

# Importing: standard
import os
import sys
import time
import tempfile
import numpy as np
from pathlib import Path
import multiprocessing
import functools
import math
import ast

# Importing: FreeCAD
import FreeCADGui as Gui
import FreeCAD as App
import Part
import Part,PartGui
import Mesh
import MeshPartGui, FreeCADGui
import MeshPart
import Mesh, Part, PartGui
import MaterialEditor
import ObjectsFem
import FemGui
import Fem
import femmesh.femmesh2mesh
from PySide import QtCore, QtGui

# Importing: paths
from freecad.chronoWorkbench                                              import ICONPATH

# Importing: util
from freecad.chronoWorkbench.util.cwloadUIfile                            import cwloadUIfile
from freecad.chronoWorkbench.util.cwloadUIicon                            import cwloadUIicon

# Importing: generation
from freecad.chronoWorkbench.generation.calc_LDPMCSL_particleVol          import calc_LDPMCSL_particleVol
from freecad.chronoWorkbench.generation.calc_LDPMCSL_surfMeshSize         import calc_LDPMCSL_surfMeshSize
from freecad.chronoWorkbench.generation.calc_LDPMCSL_surfMeshExtents      import calc_LDPMCSL_surfMeshExtents
from freecad.chronoWorkbench.generation.check_LDPMCSL_particleOverlapMPI  import check_LDPMCSL_particleOverlapMPI
from freecad.chronoWorkbench.generation.gen_CSL_facetData                 import gen_CSL_facetData
from freecad.chronoWorkbench.generation.gen_LDPMCSL_tesselation           import gen_LDPMCSL_tesselation
from freecad.chronoWorkbench.generation.gen_LDPM_facetData                import gen_LDPM_facetData
from freecad.chronoWorkbench.generation.gen_LDPMCSL_analysis              import gen_LDPMCSL_analysis
from freecad.chronoWorkbench.generation.gen_LDPMCSL_geometry              import gen_LDPMCSL_geometry
from freecad.chronoWorkbench.generation.gen_LDPMCSL_initialMesh           import gen_LDPMCSL_initialMesh
from freecad.chronoWorkbench.generation.gen_LDPMCSL_particle              import gen_LDPMCSL_particle
from freecad.chronoWorkbench.generation.gen_LDPMCSL_particleMPI           import gen_LDPMCSL_particleMPI
from freecad.chronoWorkbench.generation.gen_LDPMCSL_particleList          import gen_LDPMCSL_particleList
from freecad.chronoWorkbench.generation.gen_LDPMCSL_properties            import gen_LDPMCSL_properties
from freecad.chronoWorkbench.generation.gen_LDPMCSL_tetrahedralization    import gen_LDPMCSL_tetrahedralization

# Importing: input
from freecad.chronoWorkbench.input.read_LDPMCSL_inputs                    import read_LDPMCSL_inputs
from freecad.chronoWorkbench.input.read_LDPMCSL_tetgen                    import read_LDPMCSL_tetgen

# Importing: output
from freecad.chronoWorkbench.output.mkVtk_LDPMCSL_particles               import mkVtk_LDPMCSL_particles
from freecad.chronoWorkbench.output.mkVtk_LDPMCSL_facets                  import mkVtk_LDPMCSL_facets
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleTetFacets            import mkVtk_LDPM_singleTetFacets
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleEdgeFacets           import mkVtk_LDPM_singleEdgeFacets
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleTetParticles         import mkVtk_LDPM_singleTetParticles
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleEdgeParticles        import mkVtk_LDPM_singleEdgeParticles
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleTet                  import mkVtk_LDPM_singleTet
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleEdge                 import mkVtk_LDPM_singleEdge
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleCell                 import mkVtk_LDPM_singleCell
from freecad.chronoWorkbench.output.mkData_LDPMCSL_nodes                  import mkData_LDPMCSL_nodes
from freecad.chronoWorkbench.output.mkData_LDPMCSL_tets                   import mkData_LDPMCSL_tets
from freecad.chronoWorkbench.output.mkData_LDPMCSL_edges                  import mkData_LDPMCSL_edges
from freecad.chronoWorkbench.output.mkData_LDPMCSL_facets                 import mkData_LDPMCSL_facets
from freecad.chronoWorkbench.output.mkData_LDPMCSL_facetsVertices         import mkData_LDPMCSL_facetsVertices
from freecad.chronoWorkbench.output.mkData_LDPMCSL_particles              import mkData_LDPMCSL_particles
from freecad.chronoWorkbench.output.mkDisp_LDPMCSL_sieveCurves            import mkDisp_LDPMCSL_sieveCurves


# Turn off error for divide by zero and invalid operations
np.seterr(divide='ignore', invalid='ignore')



#sys.executable = str(Path(App.ConfigGet('AppHomePath') + '/bin/python.exe'))
multiprocessing.set_executable(str(Path(App.ConfigGet('AppHomePath') + '/bin/pythonw.exe')))



class inputWindow_LDPMCSL:
    def __init__(self):

        self.form = []

        # Load UI's for Side Panel
        self.form.append(cwloadUIfile("ui_LDPMCSL_modelProps.ui"))
        self.form.append(cwloadUIfile("ui_LDPMCSL_geometry.ui"))
        self.form.append(cwloadUIfile("ui_LDPMCSL_particles.ui"))        
        self.form.append(cwloadUIfile("ui_LDPMCSL_mixDesign.ui"))          
        self.form.append(cwloadUIfile("ui_LDPMCSL_additionalPara.ui"))       
        self.form.append(cwloadUIfile("ui_LDPMCSL_generation.ui"))

        # Label, Load Icons, and Initialize Panels
        self.form[0].setWindowTitle("Model Settings")
        self.form[1].setWindowTitle("Geometry")
        self.form[2].setWindowTitle("Particles")        
        self.form[3].setWindowTitle("Mix Design")
        self.form[4].setWindowTitle("Additional Parameters")
        self.form[5].setWindowTitle("Model Generation") 

        cwloadUIicon(self.form[0],"FEM_MaterialMechanicalNonlinear.svg")
        cwloadUIicon(self.form[1],"PartDesign_AdditiveBox.svg")
        cwloadUIicon(self.form[2],"Arch_Material_Group.svg")
        cwloadUIicon(self.form[3],"FEM_ConstraintFlowVelocity.svg")
        cwloadUIicon(self.form[4],"FEM_CreateNodesSet.svg")
        cwloadUIicon(self.form[5],"ldpm.svg")

        # Set initial output directory
        self.form[5].outputDir.setText(str(Path(App.ConfigGet('UserHomePath') + '/chronoWorkbench')))

        # Connect Open File Buttons
        QtCore.QObject.connect(self.form[0].readFileButton, QtCore.SIGNAL("clicked()"), self.openFilePara)
        QtCore.QObject.connect(self.form[1].readFileButton, QtCore.SIGNAL("clicked()"), self.openFileGeo)
        QtCore.QObject.connect(self.form[5].readDirButton, QtCore.SIGNAL("clicked()"), self.openDir)

        # Run generation for LDPM or CSL
        QtCore.QObject.connect(self.form[5].generate, QtCore.SIGNAL("clicked()"), self.generation)
        QtCore.QObject.connect(self.form[5].writePara, QtCore.SIGNAL("clicked()"), self.writeParameters)



    def getStandardButtons(self):

        # Only show a close button
        # def accept() in no longer needed, since there is no OK button
        return int(QtGui.QDialogButtonBox.Close)

    def openFilePara(self):

        path = App.ConfigGet("UserHomePath")
        filetype = "CW Parameter input format (*.cwPar)"

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getOpenFileName(None,QString.fromLocal8Bit("Read a file parameter file"),path,             filetype) # type: ignore
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getOpenFileName(None, "Read a file parameter file", path,             filetype) #PySide
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        if OpenName == "":                                                            # if the name file are not selected then Abord process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            self.form[0].setupFile.setText(OpenName)
        
        self.readParameters()

    def openFileGeo(self):

        path = App.ConfigGet("UserHomePath")
        filetype = "Supported formats (*.brep *.brp *.iges *.igs *.step *.stp);;\
                    BREP format       (*.brep *.brp);; \
                    IGES format       (*.iges *.igs);; \
                    STEP format       (*.step *.stp)"

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getOpenFileName(None,QString.fromLocal8Bit("Read a geometry file"),path,             filetype) # type: ignore
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getOpenFileName(None, "Read a geometry file", path,             filetype) #PySide
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        if OpenName == "":                                                            # if the name file are not selected then Abord process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            self.form[1].cadName.setText("Import geometry listed below")
            self.form[1].cadFile.setText(OpenName)



    def openDir(self):

        path = App.ConfigGet('UserHomePath')

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getExistingDirectory(None, "Open Directory",path,QtGui.QFileDialog.Option.ShowDirsOnly) 
         
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getExistingDirectory(None, "Open Directory",path,QtGui.QFileDialog.Option.ShowDirsOnly) 

        

        if OpenName == "":                                                            # if not selected then Abort process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            self.form[5].outputDir.setText(OpenName)

        return OpenName



    def writeParameters(self):

        # Read in inputs from input panel
        [setupFile, constitutiveEQ, matParaSet, \
            numCPU, numIncrements,maxIter,placementAlg,\
            geoType, dimensions, cadFile,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
            cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, \
            fillerC, fillerDensity, airFrac2,\
            outputDir, singleTetGen, modelType] = read_LDPMCSL_inputs(self.form)


        # Make output directory if does not exist
        outDir =  self.form[5].outputDir.text()
        try:
            os.mkdir(outDir)
        except:
            pass

        # Write parameters to file
        with open(Path(outDir + "/chronoWorkbench.cwPar"), "w") as f:
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
            f.write("outputDir = " + outputDir + "\n")

        print("Parameters written to file")

    def readParameters(self):

        paraFile = self.form[0].setupFile.text()

        # Read parameters from file
        with open(Path(paraFile), "r") as f:
            for line in f:
                if "constitutiveEQ" in line:
                    constitutiveEQ = line.split("=")[1].strip()
                elif "matParaSet" in line:
                    matParaSet = line.split("=")[1].strip()
                elif "numCPU" in line:
                    numCPU = int(line.split("=")[1].strip())
                elif "numIncrements" in line:
                    numIncrements = int(line.split("=")[1].strip())
                elif "maxIter" in line:
                    maxIter = int(line.split("=")[1].strip())
                elif "placementAlg" in line:
                    placementAlg = line.split("=")[1].strip()
                elif "geoType" in line:
                    geoType = line.split("=")[1].strip()
                elif "dimensions" in line:
                    # Change from format  ['17.00 mm', '2.00 mm', '5.00 mm'] to [17.00, 2.00, 5.00]
                    dimensions = [float(x.split()[0].strip("'")) for x in line.split("=")[1].strip().strip("[").strip("]").split(",")]
                elif "cadFile" in line:
                    cadFile = line.split("=")[1].strip()
                elif "minPar" in line:
                    minPar = float(line.split("=")[1].strip())
                elif "maxPar" in line:
                    maxPar = float(line.split("=")[1].strip())
                elif "fullerCoef" in line:
                    fullerCoef = float(line.split("=")[1].strip())
                elif "sieveCurveDiameter" in line:
                    sieveCurveDiameter = line.split("=")[1].strip()
                elif "sieveCurvePassing" in line:
                    sieveCurvePassing = line.split("=")[1].strip()
                elif "wcRatio" in line:
                    wcRatio = float(line.split("=")[1].strip())
                elif "densityWater" in line:
                    densityWater = float(line.split("=")[1].strip())
                elif "cementC" in line:
                    cementC = float(line.split("=")[1].strip())
                elif "flyashC" in line:
                    flyashC = float(line.split("=")[1].strip())
                elif "silicaC" in line:
                    silicaC = float(line.split("=")[1].strip())
                elif "scmC" in line:
                    scmC = float(line.split("=")[1].strip())
                elif "cementDensity" in line:
                    cementDensity = float(line.split("=")[1].strip())
                elif "flyashDensity" in line:
                    flyashDensity = float(line.split("=")[1].strip())
                elif "silicaDensity" in line:
                    silicaDensity = float(line.split("=")[1].strip())
                elif "scmDensity" in line:
                    scmDensity = float(line.split("=")[1].strip())
                elif "airFrac1" in line:
                    airFrac1 = float(line.split("=")[1].strip())
                elif "fillerC" in line:
                    fillerC = float(line.split("=")[1].strip())
                elif "fillerDensity" in line:
                    fillerDensity = float(line.split("=")[1].strip())
                elif "airFrac2" in line:
                    airFrac2 = float(line.split("=")[1].strip())
                elif "outputDir" in line:
                    outputDir = line.split("=")[1].strip()

        # Write parameters to input panel
        self.form[0].constEQ.setCurrentText(constitutiveEQ)
        if self.form[0].constEQ.currentIndex() == 0:
            self.form[0].matParaSet4EQ1.setCurrentText(matParaSet)
        elif self.form[0].constEQ.currentIndex() == 1:
            self.form[0].matParaSet4EQ2.setCurrentText(matParaSet)
        elif self.form[0].constEQ.currentIndex() == 2:
            self.form[0].matParaSet4EQ3.setCurrentText(matParaSet)
        elif self.form[0].constEQ.currentIndex() == 3:
            self.form[0].matParaSet4EQ4.setCurrentText(matParaSet)
        elif self.form[0].constEQ.currentIndex() == 4:
            self.form[0].matParaSet4EQ5.setCurrentText(matParaSet)
        elif self.form[0].constEQ.currentIndex() == 5:
            self.form[0].matParaSet4EQ6.setCurrentText(matParaSet)
        elif self.form[0].constEQ.currentIndex() == 6:
            self.form[0].matParaSet4EQ7.setCurrentText(matParaSet)
        elif self.form[0].constEQ.currentIndex() == 7:
            self.form[0].matParaSet4EQ8.setCurrentText(matParaSet)
        self.form[0].numCPUbox.setValue(numCPU)
        self.form[0].numPIncBox.setValue(numIncrements)
        self.form[0].numIncBox.setValue(maxIter)
        self.form[0].placementAlg.setCurrentText(placementAlg)
        self.form[1].geometryType.setCurrentText(geoType)
        if geoType == "Box":
            self.form[1].boxLength.setProperty('rawValue',(dimensions[0]))
            self.form[1].boxWidth.setProperty('rawValue',(dimensions[1]))
            self.form[1].boxHeight.setProperty('rawValue',(dimensions[2]))
        elif geoType == "Cylinder":
            self.form[1].cylinderHeight.setProperty('rawValue',(dimensions[0]))
            self.form[1].cylinderRadius.setProperty('rawValue',(dimensions[1]))
        elif geoType == "Cone":
            self.form[1].coneHeight.setProperty('rawValue',(dimensions[0]))
            self.form[1].coneRadius1.setProperty('rawValue',(dimensions[1]))
            self.form[1].coneRadius2.setProperty('rawValue',(dimensions[2]))
        elif geoType == "Sphere":
            self.form[1].sphereRadius.setProperty('rawValue',(dimensions[0]))
        elif geoType == "Ellipsoid":
            self.form[1].ellipsoidRadius1.setProperty('rawValue',(dimensions[0]))
            self.form[1].ellipsoidRadius2.setProperty('rawValue',(dimensions[1]))
            self.form[1].ellipsoidRadius3.setProperty('rawValue',(dimensions[2]))
            self.form[1].ellipsoidAngle1.setProperty('rawValue',(dimensions[3]))
            self.form[1].ellipsoidAngle2.setProperty('rawValue',(dimensions[4]))
            self.form[1].ellipsoidAngle3.setProperty('rawValue',(dimensions[5]))
        elif geoType == "Prism":
            self.form[1].prismCircumradius.setProperty('rawValue',(dimensions[0]))
            self.form[1].prismHeight.setProperty('rawValue',(dimensions[1]))
            self.form[1].prismPolygon.setProperty('rawValue',(dimensions[2]))
        elif geoType == "Notched Prism - Square":
            self.form[1].notchBoxLength.setProperty('rawValue',(dimensions[0]))
            self.form[1].notchBoxWidth.setProperty('rawValue',(dimensions[1]))
            self.form[1].notchBoxHeight.setProperty('rawValue',(dimensions[2]))
            self.form[1].notchWidth.setProperty('rawValue',(dimensions[3]))
            self.form[1].notchDepth.setProperty('rawValue',(dimensions[4]))
        elif geoType == "Notched Prism - Semi Circle":
            self.form[1].notchSCBoxLength.setProperty('rawValue',(dimensions[0]))
            self.form[1].notchSCBoxWidth.setProperty('rawValue',(dimensions[1]))
            self.form[1].notchSCBoxHeight.setProperty('rawValue',(dimensions[2]))
            self.form[1].notchSCWidth.setProperty('rawValue',(dimensions[3]))
            self.form[1].notchSCDepth.setProperty('rawValue',(dimensions[4]))
        elif geoType == "Notched Prism - Semi Ellipse":
            self.form[1].notchSEBoxLength.setProperty('rawValue',(dimensions[0]))
            self.form[1].notchSEBoxWidth.setProperty('rawValue',(dimensions[1]))
            self.form[1].notchSEBoxHeight.setProperty('rawValue',(dimensions[2]))
            self.form[1].notchSEWidth.setProperty('rawValue',(dimensions[3]))
            self.form[1].notchSEDepth.setProperty('rawValue',(dimensions[4]))
        elif geoType == "Dogbone":
            self.form[1].dogboneLength.setProperty('rawValue',(dimensions[0]))
            self.form[1].dogboneWidth.setProperty('rawValue',(dimensions[1]))
            self.form[1].dogboneThickness.setProperty('rawValue',(dimensions[2]))
            self.form[1].gaugeLength.setProperty('rawValue',(dimensions[3]))
            self.form[1].gaugeWidth.setProperty('rawValue',(dimensions[4]))
            self.form[1].dogboneType.setCurrentText(dimensions[5])
        self.form[1].cadFile.setText(cadFile)
        self.form[2].minPar.setValue(minPar)
        self.form[2].maxPar.setValue(maxPar)
        self.form[2].fullerCoef.setValue(fullerCoef)
        self.form[2].sieveDiameters.setText(str(sieveCurveDiameter))
        self.form[2].sievePassing.setText(str(sieveCurvePassing))
        self.form[3].wcRatio.setValue(wcRatio)
        self.form[3].waterDensity.setText(str(densityWater))
        self.form[3].cementContent.setText(str(cementC))
        self.form[3].flyashContent.setText(str(flyashC))
        self.form[3].silicaContent.setText(str(silicaC))
        self.form[3].scmContent.setText(str(scmC))
        self.form[3].cementDensity.setText(str(cementDensity))
        self.form[3].flyashDensity.setText(str(flyashDensity))
        self.form[3].silicaDensity.setText(str(silicaDensity))
        self.form[3].scmDensity.setText(str(scmDensity))
        self.form[3].airFrac.setValue(airFrac1)
        self.form[3].fillerContent.setText(str(fillerC))
        self.form[3].fillerDensity.setText(str(fillerDensity))
        self.form[3].airFracArb.setValue(airFrac2)
        self.form[5].outputDir.setText(outputDir)








    def generation(self):

        # Make output directory if does not exist
        outDir =  self.form[5].outputDir.text()
        try:
            os.mkdir(outDir)
        except:
            pass

        # Initialize code start time to measure performance
        start_time = time.time()

        # Make a temporary path location
        tempPath = tempfile.gettempdir() + "/chronoConc" + str(int(np.random.uniform(1e7,1e8))) + '/'
        os.mkdir(tempPath)



        # Store document
        docGui = Gui.activeDocument()

        # Make new document and set view if does not exisit
        try:
            docGui.activeView().viewAxonometric()
        except:
            App.newDocument("Unnamed")
            docGui = Gui.activeDocument()
            docGui.activeView().viewAxonometric()
        Gui.runCommand('Std_PerspectiveCamera',1)

        # Read in inputs from input panel
        [setupFile, constitutiveEQ, matParaSet, \
            numCPU, numIncrements,maxIter,placementAlg,\
            geoType, dimensions, cadFile,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
            cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, \
            fillerC, fillerDensity, airFrac2,\
            outputDir, singleTetGen, modelType] = read_LDPMCSL_inputs(self.form)

        try:
            sieveCurveDiameter = ast.literal_eval(sieveCurveDiameter)
            sieveCurvePassing = ast.literal_eval(sieveCurvePassing)
        except:
            pass

        if modelType in ["Confinement Shear Lattice (CSL) - LDPM Style ",\
                         "Confinement Shear Lattice (CSL) - Original"]:
            elementType = "CSL"
        else:
            elementType = "LDPM"


        if fillerC > 0:
            airFrac = airFrac2
        else:
            airFrac = airFrac1
        
        aggOffsetCoeff = 0.2                                    # Minimum distance between particles factor 
        verbose = "On"

        self.form[5].progressBar.setValue(1) 
        self.form[5].statusWindow.setText("Status: Generating objects.") 



        geoName = elementType + "geo" + str(0).zfill(3)
        meshName = elementType + "mesh" + str(0).zfill(3)
        analysisName = elementType + "analysis"
        materialName = elementType + "material"
        dataFilesName = elementType + 'dataFiles'+ str(0).zfill(3)
        visualFilesName = elementType + 'visualFiles'+ str(0).zfill(3)

        i = 0
        try:
            test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0] != None)
        except:
            test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName) != [])

        while test == True:
            i = i+1
            geoName = elementType + "geo" + str(i).zfill(3)
            meshName = elementType + "mesh" + str(i).zfill(3)
            dataFilesName = elementType + 'dataFiles'+ str(i).zfill(3)
            visualFilesName = elementType + 'visualFiles'+ str(i).zfill(3)
            try:
                test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0] != None)
            except:
                test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName) != [])

        # Generate geometry
        self.form[5].statusWindow.setText("Status: Generating geometry.") 
        genGeo = gen_LDPMCSL_geometry(dimensions,geoType,geoName,cadFile)
        self.form[5].progressBar.setValue(2) 

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)
        Gui.runCommand('Std_PerspectiveCamera',1)


        # Generate analysis objects
        self.form[5].statusWindow.setText("Status: Generating analysis objects.") 
        genAna = gen_LDPMCSL_analysis(analysisName,materialName)
        self.form[5].progressBar.setValue(3) 



        # Generate surface mesh
        self.form[5].statusWindow.setText("Status: Generating surface mesh.") 
        [meshVertices,meshTets, surfaceNodes,surfaceFaces] = gen_LDPMCSL_initialMesh(analysisName,geoName,meshName,minPar,maxPar)
        self.form[5].progressBar.setValue(5) 


        # Gets extents of geometry
        [minC,maxC] = calc_LDPMCSL_surfMeshExtents(meshVertices)

        self.form[5].statusWindow.setText("Status: Calculating input data.") 
        # Calculate required volume of particles and sieve curve data
        [volFracPar, tetVolume, parVolTotal,cdf,cdf1,kappa_i,newSieveCurveD, newSieveCurveP, NewSet] = calc_LDPMCSL_particleVol(wcRatio,airFrac,fullerCoef,cementC,cementDensity,densityWater,\
            flyashC,silicaC,scmC,flyashDensity,silicaDensity,scmDensity,fillerC,fillerDensity,\
            meshVertices,meshTets,minPar,maxPar,sieveCurveDiameter,sieveCurvePassing)

        # Temporary to skip over sieve curve option
        #newSieveCurveD = 0
        #NewSet = 0

        self.form[5].statusWindow.setText("Status: Calculating list of particles.") 
        # Calculate list of particle diameters for placement
        [maxParNum,parDiameterList] = gen_LDPMCSL_particleList(parVolTotal,minPar,maxPar,newSieveCurveD,\
            cdf,kappa_i,NewSet,fullerCoef)
        
        

        # Calculation of surface mesh size
        maxEdgeLength = calc_LDPMCSL_surfMeshSize(meshVertices,surfaceFaces)


        # Basic Calcs
        aggOffset = aggOffsetCoeff*minPar

        
        # Store coordinates of meshTets in new format
        coord1 = meshVertices[meshTets[:,0]-1]
        coord2 = meshVertices[meshTets[:,1]-1]
        coord3 = meshVertices[meshTets[:,2]-1]
        coord4 = meshVertices[meshTets[:,3]-1]




        verts = meshVertices[np.array(meshTets).flatten()-1]
        max_dist = np.max(np.sqrt(np.sum(verts**2, axis=1)))




        # Initialize empty particle nodes list outside geometry
        internalNodes = (np.zeros((len(parDiameterList),3))+2)*maxC




        self.form[5].statusWindow.setText('Status: Placing particles into geometry. (' + str(0) + '/' + str(len(parDiameterList)) + ')') 
        # Initialize values
        newMaxIter = 6
        particlesPlaced = 0

        




        if numCPU > 1:
        
            
            for increment in range(numIncrements-1):

                process_pool = multiprocessing.Pool(numCPU)

                outputMPI = process_pool.map(functools.partial(gen_LDPMCSL_particleMPI, surfaceNodes,maxParNum, minC, maxC, meshVertices, \
                    meshTets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                    maxPar,aggOffset,verbose,parDiameterList,maxEdgeLength,max_dist,internalNodes), parDiameterList[particlesPlaced:particlesPlaced+math.floor(len(parDiameterList)/numIncrements)])

                nodeMPI = np.array(outputMPI)[:,0:3]
                diameter = np.array(outputMPI)[:,3]
                newMaxIter = int(max(np.array(outputMPI)[:,4]))
                maxAttempts = int(max(np.array(outputMPI)[:,5]))

                particlesPlaced = particlesPlaced+len(np.array(outputMPI)[:,0:3])        

                for x in range(len(nodeMPI)):

                    # Store placed particles from this increment
                    internalNodes[particlesPlaced+x,:] = nodeMPI[x,:]

                    # Obtain extents for floating bin for node to test
                    binMin = np.array(([nodeMPI[x,0]-diameter[x]/2-maxPar/2-aggOffset,\
                        nodeMPI[x,1]-diameter[x]/2-maxPar/2-aggOffset,nodeMPI[x,2]-\
                        diameter[x]/2-maxPar/2-aggOffset]))
                    binMax = np.array(([nodeMPI[x,0]+diameter[x]/2+maxPar/2+aggOffset,\
                        nodeMPI[x,1]+diameter[x]/2+maxPar/2+aggOffset,nodeMPI[x,2]+\
                        diameter[x]/2+maxPar/2+aggOffset]))

                    # Check if particle overlapping any just added particles (ignore first one placed)
                    if x > 0:

                        overlap = check_LDPMCSL_particleOverlapMPI(nodeMPI[x,:],diameter[x],binMin,\
                            binMax,minPar,aggOffset,nodeMPI[0:x],diameter[0:x])

                        if overlap == True:

                            [newMaxIter,node,iterReq] = gen_LDPMCSL_particle(surfaceNodes,\
                                parDiameterList[particlesPlaced+x], meshVertices, \
                                meshTets,newMaxIter,maxIter,minPar,\
                                maxPar,aggOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes)
                            
                            internalNodes[particlesPlaced+x,:] = node[0,:]


                self.form[5].progressBar.setValue(95*((x)/len(parDiameterList))+6) 
                self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')


        # Generate particles for length of needed aggregate (not placed via MPI)
        for x in range(particlesPlaced,len(parDiameterList)):

            # Generate particle
            [newMaxIter,node,iterReq] = gen_LDPMCSL_particle(surfaceNodes,parDiameterList[x],meshVertices,meshTets,newMaxIter,maxIter,minPar,maxPar,\
                aggOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes)

            # Update progress bar every 1% of placement
            if x % np.rint(len(parDiameterList)/100) == 0:
                self.form[5].progressBar.setValue(80*((x)/len(parDiameterList))+6) 

            if len(parDiameterList)<=1000:
                # Update number particles placed every 1%
                if x % np.rint(len(parDiameterList)/100) == 0:
                    self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')
            elif len(parDiameterList)<=10000:
                # Update number particles placed every 0.1%
                if x % np.rint(len(parDiameterList)/1000) == 0:
                    self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')
            else:
                # Update number particles placed every 0.01%
                if x % np.rint(len(parDiameterList)/10000) == 0:
                    self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')

            internalNodes[x,:] = node

        self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(len(parDiameterList)) + '/' + str(len(parDiameterList)) + ')')


        materialList = np.ones(len(parDiameterList))

        placementTime = round(time.time() - start_time,2)   
        nParticles = len(parDiameterList)

        # Create empty lists if not multi-material or cementStructure
        aggGrainsDiameterList, itzDiameterList, binderDiameterList, PoresDiameterList,\
            ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0,0,0,0









        tetTessTimeStart = time.time()

        # Generate tetrahedralization
        self.form[5].statusWindow.setText("Status: Forming tetrahedralization.") 
        tetGen = gen_LDPMCSL_tetrahedralization(internalNodes,surfaceNodes,\
            surfaceFaces,geoName,tempPath)
        self.form[5].progressBar.setValue(89) 


        # Read in tetrahedralization
        [allNodes,allTets,allEdges] = read_LDPMCSL_tetgen(Path(tempPath + geoName \
        + '.node'),Path(tempPath + geoName + '.ele'),Path(tempPath + geoName + '.edge'))
        self.form[5].progressBar.setValue(90) 


        # Generate tesselation
        self.form[5].statusWindow.setText("Status: Forming tesselation.") 
        

        [tetFacets,facetCenters,facetAreas,facetNormals,tetn1,tetn2,tetPoints,allDiameters,facetPointData,facetCellData] = \
            gen_LDPMCSL_tesselation(allNodes,allTets,parDiameterList,minPar,geoName)    






        self.form[5].progressBar.setValue(95) 
        tetTessTime = round(time.time() - tetTessTimeStart,2)   




        # Store values for unused features
        edgeMaterialList = 0
        materialRule = 0
        multiMaterial = 'Off'
        cementStructure = 'Off'




        writeTimeStart = time.time()





     
        self.form[5].statusWindow.setText("Status: Generating facet data information.") 

        if elementType == "LDPM":
            [facetData,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial] = gen_LDPM_facetData(\
                allNodes,allTets,tetFacets,facetCenters,facetAreas,facetNormals,tetn1,\
                tetn2,materialList,materialRule,multiMaterial,cementStructure,edgeMaterialList,facetCellData)
        elif elementType == "CSL":
            [facetData,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial] = gen_CSL_facetData(\
                allNodes,allEdges,allTets,tetFacets,facetCenters,facetAreas,facetNormals,tetn1,\
                tetn2,materialList,materialRule,multiMaterial,cementStructure,edgeMaterialList,facetCellData)


        self.form[5].progressBar.setValue(98) 







        self.form[5].statusWindow.setText("Status: Writing external facet data file.") 
        # Create file of external triangle facets for plotting of cells
        #externalFacetsFile = externalFacetFile(facetData,meshVertices,surfaceFaces,geoName)





        # Initialize counter for number of facet materials switched
        matSwitched = 0

        itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct,\
            PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
            PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,\
            matSwitched,materialRule = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0




        App.activeDocument().addObject('App::DocumentObjectGroup',dataFilesName)
        App.activeDocument().getObject(dataFilesName).Label = 'Data Files'

        App.activeDocument().addObject('App::DocumentObjectGroup',visualFilesName)
        App.activeDocument().getObject(visualFilesName).Label = 'Visualization Files'




        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName))
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName))


        self.form[5].statusWindow.setText("Status: Writing node data file.")

        mkData_LDPMCSL_nodes(geoName,tempPath,allNodes)

        self.form[5].statusWindow.setText("Status: Writing tet data file.")

        mkData_LDPMCSL_tets(geoName,tempPath,allTets)

        if elementType == "CSL":
            self.form[5].statusWindow.setText("Status: Writing edge data file.")
            mkData_LDPMCSL_edges(geoName,tempPath,allEdges)


        self.form[5].statusWindow.setText("Status: Writing facet data file.")

        # If data files requested, generate Facet File
        mkData_LDPMCSL_facets(geoName,tempPath,facetData)
        mkData_LDPMCSL_facetsVertices(geoName,tempPath,tetFacets)


        self.form[5].statusWindow.setText("Status: Writing particle data file.")

        # If data files requested, generate Particle Data File
        mkData_LDPMCSL_particles(allNodes,parDiameterList,geoName,tempPath)


        self.form[5].statusWindow.setText("Status: Writing visualization files.")





        # If visuals requested, generate Particle VTK File
        mkVtk_LDPMCSL_particles(internalNodes,parDiameterList,materialList,geoName,tempPath)

        # If visuals requested, generate Facet VTK File
        mkVtk_LDPMCSL_facets(geoName,tempPath,facetPointData,facetCellData)
        
        if singleTetGen == True:
            if elementType == "LDPM":
                mkVtk_LDPM_singleTetFacets(geoName,tempPath,tetFacets)
                mkVtk_LDPM_singleTetParticles(allNodes,allTets,allDiameters,geoName,tempPath)
                mkVtk_LDPM_singleTet(allNodes,allTets,geoName,tempPath)
                mkVtk_LDPM_singleCell(allNodes,allTets,parDiameterList,tetFacets,geoName,tempPath)
            elif elementType == "CSL":
                pass
                mkVtk_LDPM_singleEdgeFacets(geoName,tempPath,allEdges,facetData,tetFacets)
                mkVtk_LDPM_singleEdgeParticles(allNodes,allEdges,allDiameters,geoName,tempPath)
                mkVtk_LDPM_singleEdge(allNodes,allEdges,geoName,tempPath)
            else:
                pass

        writeTime = round(time.time() - writeTimeStart,2)



        # Generate Log file after run
        #mkLogFile = logFile(gmshTime,nParticles,placementTime,maxPar,\
        #    minPar,fullerCoef,wcRatio,cementC,volFracAir,q,maxIter,\
        #    geoName,aggOffset,densityWater,densityCement,allTets,dataType,\
        #    tetTessTime,writeTime,geoFile,dFiber,lFiber,vFiber,fiberFile,\
        #    multiMaterial,materialFile,maxGrainD,minGrainD,grainFullerCoef,\
        #    maxBinderD,minBinderD,binderFullerCoef,maxITZD,minITZD,ITZFullerCoef,output,fibers,\
        #    itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
        #    aggVolFracAct,sieveCurveDiameter,sieveCurvePassing,matSwitched,materialRule,\
        #    cementStructure,cementmaterialFile,maxPoresD,minPoresD,PoresFullerCoef,\
        #    PoresSieveCurveDiameter,PoresSieveCurvePassing,maxClinkerD,minClinkerD,\
        #    ClinkerFullerCoef,ClinkerSieveCurveDiameter,ClinkerSieveCurvePassing,\
        #    maxCHD,minCHD,CHFullerCoef,CHSieveCurveDiameter,CHSieveCurvePassing,\
        #    maxCSH_LDD,minCSH_LDD,CSH_LDFullerCoef,CSH_LDSieveCurveDiameter,CSH_LDSieveCurvePassing,\
        #    maxCSH_HDD,minCSH_HDD,CSH_HDFullerCoef,CSH_HDSieveCurveDiameter,CSH_HDSieveCurvePassing,\
        #    PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
        #    PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,outputUnits)





        # Move files to selected output directory
        print('Moving files.')


        outName = '/' + geoName + geoType + str(i).zfill(3)
        i = 0
        while os.path.isdir(Path(outDir + outName)):
            i = i+1
            outName = '/' + geoName + geoType + str(i).zfill(3)
            
        os.rename(Path(tempPath),Path(outDir + outName))
        os.rename(Path(outDir + outName + '/' + geoName + '-para-mesh.vtk'),Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))
        os.remove(Path(outDir + outName + '/' + geoName + '2D.mesh'))
        os.remove(Path(outDir + outName + '/' + geoName + '.node'))
        os.remove(Path(outDir + outName + '/' + geoName + '.ele'))
        os.remove(Path(outDir + outName + '/' + geoName + '.edge'))



        print("Generated files written to: " + str(Path(outDir + outName)))





        # Set linked object for node data file
        LDPMnodesData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMnodesData")                                     # create your object
        LDPMnodesData.ViewObject.Proxy = IconViewProviderToFile(LDPMnodesData,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName).addObject(LDPMnodesData)
        LDPMnodesData.addProperty("App::PropertyFile",'Location','Node Data File','Location of node data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-nodes.dat'))
        
        # Set linked object for tet data file
        LDPMtetsData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMtetsData")                                     # create your object
        LDPMtetsData.ViewObject.Proxy = IconViewProviderToFile(LDPMtetsData,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName).addObject(LDPMtetsData)
        LDPMtetsData.addProperty("App::PropertyFile",'Location','Tet Data File','Location of tets data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-tets.dat'))

        # Set linked object for facet data file
        LDPMfacetsData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMfacetsData")                                     # create your object
        LDPMfacetsData.ViewObject.Proxy = IconViewProviderToFile(LDPMfacetsData,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName).addObject(LDPMfacetsData)
        LDPMfacetsData.addProperty("App::PropertyFile",'Location','Facet Data File','Location of facet data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-facets.dat'))


        # Set linked object for particle VTK file
        LDPMparticlesVTK = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMparticlesVTK")                                     # create your object
        LDPMparticlesVTK.ViewObject.Proxy = IconViewProviderToFile(LDPMparticlesVTK,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMparticlesVTK)
        LDPMparticlesVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-particles.000.vtk'))


        if geoType in ['Cylinder', 'Cone', 'Sphere']:


            # Insert mesh visualization and link mesh VTK file
            Fem.insert(str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk')),App.ActiveDocument.Name)        
            LDPMmeshVTK = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000')
            LDPMmeshVTK.Label = 'LDPMmeshVTK' 
            App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMmeshVTK)
            LDPMmeshVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))

            Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').ShapeColor = (0.80,0.80,0.80)
            Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').BackfaceCulling = False
            Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').MaxFacesShowInner = 0
            Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').DisplayMode = u"Faces"

            objName = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0].Name
            try:
                Gui.getDocument(App.ActiveDocument.Name).getObject(objName).Transparency = 50
                Gui.getDocument(App.ActiveDocument.Name).getObject(objName).ShapeColor = (0.80,0.80,0.80)
            except:
                pass

        else:
            # Set linked object for mesh VTK file
            LDPMmeshVTK = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMmeshVTK")                                     # create your object
            LDPMmeshVTK.ViewObject.Proxy = IconViewProviderToFile(LDPMmeshVTK,os.path.join(ICONPATH,'FEMMeshICON.svg'))
            App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMmeshVTK)
            LDPMmeshVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))

            objName = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0].Name
            try:
                Gui.getDocument(App.ActiveDocument.Name).getObject(objName).Transparency = 0
                Gui.getDocument(App.ActiveDocument.Name).getObject(objName).ShapeColor = (0.80,0.80,0.80)
            except:
                pass


        # Insert facet visualization and link facet VTK file
        Fem.insert(str(Path(outDir + outName + '/' + geoName + '-para-facets.000.vtk')),App.ActiveDocument.Name)        
        LDPMfacetsVTK = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000')
        LDPMfacetsVTK.Label = 'LDPMfacetsVTK' 
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMfacetsVTK)
        LDPMfacetsVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-facets.000.vtk'))


        # Set visualization properties for facets
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').DisplayMode = u"Wireframe"
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').MaxFacesShowInner = 0
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').BackfaceCulling = False
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').ShapeColor = (0.36,0.36,0.36)


        # Set visualization properties for particle centers
        Gui.getDocument(App.ActiveDocument.Name).getObject(meshName).DisplayMode = u"Nodes"
        Gui.getDocument(App.ActiveDocument.Name).getObject(meshName).PointSize = 3.00
        Gui.getDocument(App.ActiveDocument.Name).getObject(meshName).PointColor = (0.00,0.00,0.00)









        # Set material properties, descriptions, and values based on selected constitutive equation set and material parameter set
        [materialProps, materialPropDesc, materialPropsVal] = gen_LDPMCSL_properties(constitutiveEQ, matParaSet)




        simProps = [\
            "TotalTime",\
            "TimestepSize",\
            "NumberOfThreads",\
            "NumberOutputSteps",\
            ]




        simPropDesc = [\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            ]





        # Remove unused material properties
        #App.getDocument(App.ActiveDocument.Name).getObject(materialName).removeProperty("References")

        # Add appropriate constitutive equation set
        App.getDocument(App.ActiveDocument.Name).getObject(materialName).addProperty("App::PropertyString",'ConstitutiveEquationSet','Base','Set of constitutive equations.').ConstitutiveEquationSet=constitutiveEQ



        # Add appropriate material properties
        for x in range(len(materialProps)):
            App.getDocument(App.ActiveDocument.Name).getObject(materialName).addProperty("App::PropertyString",materialProps[x],elementType+" Parameters",materialPropDesc[x])
            setattr(App.getDocument(App.ActiveDocument.Name).getObject(materialName),materialProps[x],str(materialPropsVal[x]))


        # Add appropriate simulation properties
        for x in range(len(simProps)):
            App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyFloat",simProps[x],"Simulation",simPropDesc[x])#.Density=0.25
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyEnumeration","Solver","Simulation","Solver software").Solver=['Project Chrono']
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyEnumeration","IntegrationScheme","Simulation","Integrator type").IntegrationScheme=['Explicit']












        self.form[5].progressBar.setValue(100) 



        # Display sieve curve data
        mkDisp_LDPMCSL_sieveCurves(volFracPar, tetVolume, minPar, maxPar,fullerCoef,sieveCurveDiameter,sieveCurvePassing,parDiameterList)

        # Switch back to model window
        mw=Gui.getMainWindow()
        mdi=mw.findChild(QtGui.QMdiArea)
        mdi.activatePreviousSubWindow()

        # Switch to FEM GUI
        App.ActiveDocument.recompute()



        Gui.Control.closeDialog()
        Gui.activateWorkbench("FemWorkbench")
        FemGui.setActiveAnalysis(App.activeDocument().getObject(analysisName))
        







        
    # What to do when "Close" Button Clicked
    def reject(self):
         try:
             Gui.ActiveDocument.resetEdit()
             Gui.Control.closeDialog()
         except:
             Gui.Control.closeDialog()

        
class IconViewProviderToFile:                                       # Class ViewProvider create Property view of object
    def __init__( self, obj, icon):
        self.icone = icon
        
    def getIcon(self):                                              # GetIcon
        return self.icone        


class input_LDPMCSL_Class():
    """My new command"""

    def GetResources(self):
        return {"Pixmap"  : os.path.join(ICONPATH, "ldpm.svg"), # the name of a svg file available in the resources
                "MenuText": "LDPM/CSL Generation",
                "ToolTip" : "Generation of an LDPM or CSL geometry"}

    def Activated(self):

        Gui.Control.showDialog(inputWindow_LDPMCSL())

        return

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True

Gui.addCommand("mod_LDPMCSL", input_LDPMCSL_Class())