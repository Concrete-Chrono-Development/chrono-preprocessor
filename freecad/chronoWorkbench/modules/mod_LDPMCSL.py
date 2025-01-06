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
from freecad.chronoWorkbench.generation.driver_LDPMCSL                    import driver_LDPMCSL
from freecad.chronoWorkbench.generation.gen_LDPMCSL_tesselation           import gen_LDPMCSL_tesselation
from freecad.chronoWorkbench.generation.gen_LDPM_facetData                import gen_LDPM_facetData
from freecad.chronoWorkbench.generation.gen_LDPM_debugTet                 import gen_LDPM_debugTet
from freecad.chronoWorkbench.generation.gen_LDPMCSL_analysis              import gen_LDPMCSL_analysis


# Importing: input
from freecad.chronoWorkbench.input.read_LDPMCSL_inputs                    import read_LDPMCSL_inputs

# Importing: output
from freecad.chronoWorkbench.output.mkVtk_particles                       import mkVtk_particles
from freecad.chronoWorkbench.output.mkVtk_LDPMCSL_facets                  import mkVtk_LDPMCSL_facets
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleTet                  import mkVtk_LDPM_singleTet
from freecad.chronoWorkbench.output.mkPy_LDPM_singleDebugParaview         import mkPy_LDPM_singleDebugParaview
from freecad.chronoWorkbench.output.mkPy_LDPM_singleDebugParaviewLabels   import mkPy_LDPM_singleDebugParaviewLabels
from freecad.chronoWorkbench.output.mkData_nodes                          import mkData_nodes
from freecad.chronoWorkbench.output.mkData_LDPMCSL_tets                   import mkData_LDPMCSL_tets
from freecad.chronoWorkbench.output.mkData_LDPMCSL_facets                 import mkData_LDPMCSL_facets
from freecad.chronoWorkbench.output.mkData_LDPMCSL_facetsVertices         import mkData_LDPMCSL_facetsVertices
from freecad.chronoWorkbench.output.mkData_particles                      import mkData_particles
from freecad.chronoWorkbench.output.mkParameters                          import mkParameters



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
        self.form.append(cwloadUIfile("ui_LDPMCSL_debugging.ui"))

        # Label, Load Icons, and Initialize Panels
        self.form[0].setWindowTitle("Model Settings")
        self.form[1].setWindowTitle("Geometry")
        self.form[2].setWindowTitle("Particles")        
        self.form[3].setWindowTitle("Mix Design")
        self.form[4].setWindowTitle("Additional Parameters")
        self.form[5].setWindowTitle("Model Generation") 
        self.form[6].setWindowTitle("Debugging")

        cwloadUIicon(self.form[0],"FEM_MaterialMechanicalNonlinear.svg")
        cwloadUIicon(self.form[1],"PartDesign_AdditiveBox.svg")
        cwloadUIicon(self.form[2],"Arch_Material_Group.svg")
        cwloadUIicon(self.form[3],"FEM_ConstraintFlowVelocity.svg")
        cwloadUIicon(self.form[4],"FEM_CreateNodesSet.svg")
        cwloadUIicon(self.form[5],"ldpm.svg")
        cwloadUIicon(self.form[6],"PartDesign_Sprocket.svg")

        # Set initial output directory
        self.form[5].outputDir.setText(str(Path(App.ConfigGet('UserHomePath') + '/chronoWorkbench')))

        # Connect selectObject and deselectObject Button in the Geometry Tab
        QtCore.QObject.connect(self.form[1].selectObjectButton, QtCore.SIGNAL("clicked()"), self.selectGeometry)
        QtCore.QObject.connect(self.form[1].deselectObjectButton, QtCore.SIGNAL("clicked()"), self.deselectGeometry)

        # Connect Open File Buttons
        QtCore.QObject.connect(self.form[0].readFileButton, QtCore.SIGNAL("clicked()"), self.openFilePara)
        QtCore.QObject.connect(self.form[1].readFileButton, QtCore.SIGNAL("clicked()"), self.openFileGeo)
        QtCore.QObject.connect(self.form[5].readDirButton, QtCore.SIGNAL("clicked()"), self.openDir)
        QtCore.QObject.connect(self.form[4].readMultiMatFile, QtCore.SIGNAL("clicked()"), self.openMultiMatFile)
        QtCore.QObject.connect(self.form[4].readAggFile, QtCore.SIGNAL("clicked()"), self.openAggFile)

        # Run generation for LDPM or CSL
        QtCore.QObject.connect(self.form[5].generate, QtCore.SIGNAL("clicked()"), self.generationDriver)
        QtCore.QObject.connect(self.form[5].generateFast, QtCore.SIGNAL("clicked()"), self.generationDriverFast)
        QtCore.QObject.connect(self.form[5].writePara, QtCore.SIGNAL("clicked()"), self.writeParameters)

        # Run debugging generation of single tetrahedron
        QtCore.QObject.connect(self.form[6].generate_reg, QtCore.SIGNAL("clicked()"), self.debugGenerateRegTet)
        QtCore.QObject.connect(self.form[6].generate_irreg, QtCore.SIGNAL("clicked()"), self.debugGenerateIrregTet)



    def generationDriver(self):

        # Make a temporary path location
        tempPath = tempfile.gettempdir() + "/chronoConc" + str(int(np.random.uniform(1e7,1e8))) + '/'
        os.mkdir(tempPath)

        fastGen = False
        mkParameters(self,"LDPMCSL",tempPath)
        driver_LDPMCSL(self,fastGen,tempPath)

    def generationDriverFast(self):

        # Make a temporary path location
        tempPath = tempfile.gettempdir() + "/chronoConc" + str(int(np.random.uniform(1e7,1e8))) + '/'
        os.mkdir(tempPath)

        fastGen = True
        mkParameters(self,"LDPMCSL",tempPath)
        driver_LDPMCSL(self,fastGen,tempPath)


    def getStandardButtons(self):

        # Only show a close button
        # def accept() in no longer needed, since there is no OK button
        return int(QtGui.QDialogButtonBox.Close)

    def selectGeometry(self):
            
        # Select Geometry
        sel = FreeCADGui.Selection.getSelection()
        if len(sel) == 0:
            App.Console.PrintMessage("Please select a geometry"+"\n")
        else:
            self.form[1].selectedObject.setText(sel[0].Name)

    def deselectGeometry(self):
                
        # Deselect Geometry
        self.form[1].selectedObject.setText("")


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
        filetype = "Supported formats (*.brep *.brp *.iges *.igs *.step *.stp *.inp *.vtk *.vtu);;\
                    BREP format       (*.brep *.brp);; \
                    IGES format       (*.iges *.igs);; \
                    STEP format       (*.step *.stp);; \
                    Abaqus/CalcuLix format  (*.inp);; \
                    VTK Legacy/modern format (*.vtk *.vtu)"

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


    def openMultiMatFile(self):

        path = App.ConfigGet("UserHomePath")
        filetype = "Voxel Data File (*.img)"

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getOpenFileName(None,QString.fromLocal8Bit("Read a voxel data file"),path,             filetype) # type: ignore
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getOpenFileName(None, "Read a voxel data file", path,             filetype) #PySide
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        if OpenName == "":                                                            # if the name file are not selected then Abord process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            self.form[4].multiMatFile.setText(OpenName)

    def openAggFile(self):

        path = App.ConfigGet("UserHomePath")
        filetype = "Voxel Data File (*.pimg)"

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getOpenFileName(None,QString.fromLocal8Bit("Read a voxel data file"),path,             filetype) # type: ignore
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getOpenFileName(None, "Read a voxel data file", path,             filetype) #PySide
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        if OpenName == "":                                                            # if the name file are not selected then Abord process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            self.form[4].aggFile.setText(OpenName)

    def writeParameters(self):

        mkParameters(self,"LDPMCSL","writeOnly")

    # Read parameters from file and write to input panel in FreeCAD
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
                    if line != []:
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
                elif "htcToggle" in line:
                    htcToggle = line.split("=")[1].strip()
                elif "htcLength" in line:
                    htcLength = float(line.split("=")[1].strip())
                elif "fiberToggle" in line:
                    fiberToggle = line.split("=")[1].strip()
                elif "fiberCutting" in line:
                    fiberCutting = line.split("=")[1].strip()
                elif "fiberDiameter" in line:
                    fiberDiameter = float(line.split("=")[1].strip())
                elif "fiberLength" in line:
                    fiberLength = float(line.split("=")[1].strip())
                elif "fiberVol" in line:
                    fiberVol = float(line.split("=")[1].strip())
                elif "fiberOrientation1" in line:
                    fiberOrientation1 = float(line.split("=")[1].strip())
                elif "fiberOrientation2" in line:
                    fiberOrientation2 = float(line.split("=")[1].strip())
                elif "fiberOrientation3" in line:
                    fiberOrientation3 = float(line.split("=")[1].strip())
                elif "fiberPref" in line:
                    fiberPref = float(line.split("=")[1].strip())
                elif "fiberFile" in line:
                    fiberFile = line.split("=")[1].strip()
                elif "multiMatToggle" in line:
                    multiMatToggle = line.split("=")[1].strip()
                elif "aggFile" in line:
                    aggFile = line.split("=")[1].strip()
                elif "multiMatFile" in line:
                    multiMatFile = line.split("=")[1].strip()
                elif "multiMatRule" in line:
                    multiMatRule = line.split("=")[1].strip()
                elif "grainAggMin" in line:
                    grainAggMin = float(line.split("=")[1].strip())
                elif "grainAggMax" in line:
                    grainAggMax = float(line.split("=")[1].strip())
                elif "grainAggFuller" in line:
                    grainAggFuller = float(line.split("=")[1].strip())
                elif "grainAggSieveD" in line:
                    grainAggSieveD = line.split("=")[1].strip()
                elif "grainAggSieveP" in line:
                    grainAggSieveP = line.split("=")[1].strip()
                elif "grainITZMin" in line:
                    grainITZMin = float(line.split("=")[1].strip())
                elif "grainITZMax" in line:
                    grainITZMax = float(line.split("=")[1].strip())
                elif "grainITZFuller" in line:
                    grainITZFuller = float(line.split("=")[1].strip())
                elif "grainITZSieveD" in line:
                    grainITZSieveD = line.split("=")[1].strip()
                elif "grainITZSieveP" in line:
                    grainITZSieveP = line.split("=")[1].strip()
                elif "grainBinderMin" in line:
                    grainBinderMin = float(line.split("=")[1].strip())
                elif "grainBinderMax" in line:
                    grainBinderMax = float(line.split("=")[1].strip())
                elif "grainBinderFuller" in line:
                    grainBinderFuller = float(line.split("=")[1].strip())
                elif "grainBinderSieveD" in line:
                    grainBinderSieveD = line.split("=")[1].strip()
                elif "grainBinderSieveP" in line:
                    grainBinderSieveP = line.split("=")[1].strip()
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
        elif geoType == "Arbitrary Prism":
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
        if htcToggle == "On":
            self.form[4].PrimitiveTypeCB.setCurrentText("HTC Parameters")
            self.form[4].widgetStack2.setCurrentIndex(0) # Index 0 is the current page for HTC
        self.form[4].HTCtoggle.setCurrentText(htcToggle)
        self.form[4].HTClength.setValue(htcLength)
        if fiberToggle == "On":
            self.form[4].PrimitiveTypeCB.setCurrentText("Fiber-Reinforcement Parameters")
            self.form[4].widgetStack2.setCurrentIndex(1) # Index 1 is the current page for Fiber Parameters
        self.form[4].fiberToggle.setCurrentText(fiberToggle)
        self.form[4].fiberCutting.setCurrentText(fiberCutting)
        self.form[4].fiberDiameter.setValue(fiberDiameter)
        self.form[4].fiberLength.setValue(fiberLength)
        self.form[4].fiberVol.setValue(fiberVol)
        self.form[4].fiberOrien1.setValue(fiberOrientation1)
        self.form[4].fiberOrien2.setValue(fiberOrientation2)
        self.form[4].fiberOrien3.setValue(fiberOrientation3)
        self.form[4].fiberPref.setValue(fiberPref)
        self.form[4].fiberFile.setText(fiberFile)           
        if multiMatToggle == "On":
            self.form[4].PrimitiveTypeCB.setCurrentText("Multi-Material Parameters")
            self.form[4].widgetStack2.setCurrentIndex(4) # Index 4 is the current page for multi-material (P-LDPM setup)  
        self.form[4].multiMatToggle.setCurrentText(multiMatToggle)
        self.form[4].aggFile.setText(aggFile)
        self.form[4].multiMatFile.setText(multiMatFile)
        self.form[4].multiMatRule.setValue(int(multiMatRule))
        self.form[4].grainAggMin.setValue(grainAggMin)
        self.form[4].grainAggMax.setValue(grainAggMax)
        self.form[4].grainAggFuller.setValue(grainAggFuller)
        self.form[4].grainAggSieveD.setText(str(grainAggSieveD))
        self.form[4].grainAggSieveP.setText(str(grainAggSieveP))
        self.form[4].grainITZMin.setValue(grainITZMin)
        self.form[4].grainITZMax.setValue(grainITZMax)
        self.form[4].grainITZFuller.setValue(grainITZFuller)
        self.form[4].grainITZSieveD.setText(str(grainITZSieveD))
        self.form[4].grainITZSieveP.setText(str(grainITZSieveP))
        self.form[4].grainBinderMin.setValue(grainBinderMin)
        self.form[4].grainBinderMax.setValue(grainBinderMax)
        self.form[4].grainBinderFuller.setValue(grainBinderFuller)
        self.form[4].grainBinderSieveD.setText(str(grainBinderSieveD))
        self.form[4].grainBinderSieveP.setText(str(grainBinderSieveP))
        self.form[5].outputDir.setText(outputDir)


    def debugGenerateRegTet(self):

        self.debugGenerateTet("Regular")

    def debugGenerateIrregTet(self):

        self.debugGenerateTet("Irregular")


    def debugGenerateTet(self,type):


        # Make output directory if does not exist
        outDir =  self.form[5].outputDir.text()
        try:
            os.mkdir(outDir)
        except:
            pass

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

        # Read in inputs from input panel just to get location for file output
        [setupFile, constitutiveEQ, matParaSet, \
            numCPU, numIncrements,maxIter,placementAlg,\
            geoType, dimensions, cadFile,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
            cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, \
            fillerC, fillerDensity, airFrac2,\
            htcToggle, htcLength,\
            fiberToggle, fiberCutting, fiberDiameter, fiberLength, fiberVol, fiberOrientation1, fiberOrientation2, fiberOrientation3, fiberPref, fiberFile,\
            multiMatToggle,aggFile,multiMatFile,multiMatRule,\
            grainAggMin, grainAggMax, grainAggFuller, grainAggSieveD, grainAggSieveP,\
            grainITZMin, grainITZMax, grainITZFuller, grainITZSieveD, grainITZSieveP,\
            grainBinderMin, grainBinderMax, grainBinderFuller, grainBinderSieveD, grainBinderSieveP,\
            periodicToggle,\
            outputDir, dataFilesGen, visFilesGen, singleTetGen, modelType] = read_LDPMCSL_inputs(self.form)

        if modelType in ["Confinement Shear Lattice (CSL) - LDPM Style ",\
                         "Confinement Shear Lattice (CSL) - Original"]:
            elementType = "CSL"
        else:
            elementType = "LDPM"

        geoName = elementType + "geo" + str(0).zfill(3)
        meshName = elementType + "mesh" + str(0).zfill(3)
        analysisName = elementType + "analysis"
        materialName = elementType + "material"
        dataFilesName = elementType + 'dataFiles'+ str(0).zfill(3)
        visualFilesName = elementType + 'visualFiles'+ str(0).zfill(3)

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)
        Gui.runCommand('Std_PerspectiveCamera',1)


        # Generate analysis objects
        self.form[5].statusWindow.setText("Status: Generating analysis objects.") 
        genAna = gen_LDPMCSL_analysis(analysisName,materialName)
        self.form[5].progressBar.setValue(3) 

        [allNodes,allTets,parDiameterList,materialList,minPar,geoName] = gen_LDPM_debugTet(type)

        [tetFacets,facetCenters,facetAreas,facetNormals,tetn1,tetn2,tetPoints,allDiameters,facetPointData,facetCellData] = \
            gen_LDPMCSL_tesselation(allNodes,allTets,parDiameterList,minPar,geoName)   


        self.form[5].progressBar.setValue(95) 

        # Store values for unused features
        edgedata = 0
        edgeMaterialList = 0
        multiMatRule = 0
        particleID = np.zeros(4)
        multiMaterial = 'Off'
        cementStructure = 'Off'

        [facetData,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial] = gen_LDPM_facetData(\
            allNodes,allTets,tetFacets,facetCenters,facetAreas,facetNormals,tetn1,\
            tetn2,materialList,multiMatRule,multiMaterial,cementStructure,edgeMaterialList,facetCellData,particleID)

        self.form[5].progressBar.setValue(98) 


        self.form[5].statusWindow.setText("Status: Writing external facet data file.") 

        # Initialize counter for number of facet materials switched
        matSwitched = 0

        itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct,\
            PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
            PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,\
            matSwitched,multiMatRule = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0


        App.activeDocument().addObject('App::DocumentObjectGroup',dataFilesName)
        App.activeDocument().getObject(dataFilesName).Label = 'Data Files'

        App.activeDocument().addObject('App::DocumentObjectGroup',visualFilesName)
        App.activeDocument().getObject(visualFilesName).Label = 'Visualization Files'

        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName))
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName))

        self.form[5].statusWindow.setText("Status: Writing node data file.")

        mkData_nodes(geoName,tempPath,allNodes)

        self.form[5].statusWindow.setText("Status: Writing tet data file.")

        mkData_LDPMCSL_tets(geoName,tempPath,allTets)

        self.form[5].statusWindow.setText("Status: Writing facet data file.")

        # If data files requested, generate Facet File
        mkData_LDPMCSL_facets(geoName,tempPath,facetData)
        mkData_LDPMCSL_facetsVertices(geoName,tempPath,tetFacets)
        #mkData_LDPMCSL_faceFacets(geoName,tempPath,surfaceNodes,surfaceFaces)

        self.form[5].statusWindow.setText("Status: Writing particle data file.")

        # If data files requested, generate Particle Data File
        mkData_particles(allNodes,parDiameterList,geoName,tempPath)

        self.form[5].statusWindow.setText("Status: Writing visualization files.")

        # If visuals requested, generate Particle VTK File
        mkVtk_particles(allNodes,parDiameterList,materialList,geoName,tempPath)

        # If visuals requested, generate Facet VTK File
        mkVtk_LDPMCSL_facets(geoName,tempPath,tetFacets,facetMaterial)

        mkVtk_LDPM_singleTet(allNodes,allTets,geoName,tempPath)


        i = 0
        outName = '/' + geoName + geoType + str(i).zfill(3)
        while os.path.isdir(Path(outDir + outName)):
            i = i+1
            outName = '/' + geoName + geoType + str(i).zfill(3)


        mkPy_LDPM_singleDebugParaview(geoName, outDir, outName, tempPath)
        mkPy_LDPM_singleDebugParaviewLabels(geoName, tempPath)

        # Move files to selected output directory
        print('Moving files.')

           
        os.rename(Path(tempPath),Path(outDir + outName))
        os.rename(Path(outDir + outName + '/' + geoName + '-para-singleTet.000.vtk'),Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))



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



        #objName = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0].Name
        #try:
        #    Gui.getDocument(App.ActiveDocument.Name).getObject(objName).Transparency = 0
        #    Gui.getDocument(App.ActiveDocument.Name).getObject(objName).ShapeColor = (0.80,0.80,0.80)
        #except:
        #    pass



        # Insert mesh visualization and link mesh VTK file
        Fem.insert(str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk')),App.ActiveDocument.Name)
        LDPMmeshVTK = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000')
        LDPMmeshVTK.Label = 'LDPMmeshVTK'
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMmeshVTK)
        LDPMmeshVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))

        # Set visualization properties for mesh
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').DisplayMode = u"Wireframe & Nodes"
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').PointSize = 20.00
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').PointColor = (255.0,170.0,0.0)
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').MaxFacesShowInner = 0
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').BackfaceCulling = False
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').ShapeColor = (0.36,0.36,0.36)
        









        # Insert facet visualization and link facet VTK file
        Fem.insert(str(Path(outDir + outName + '/' + geoName + '-para-facets.000.vtk')),App.ActiveDocument.Name)        
        LDPMfacetsVTK = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000')
        LDPMfacetsVTK.Label = 'LDPMfacetsVTK' 
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMfacetsVTK)
        LDPMfacetsVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-facets.000.vtk'))


        # Set visualization properties for facets
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').DisplayMode = u"Faces & Wireframe"
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').MaxFacesShowInner = 0
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').BackfaceCulling = False
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').ShapeColor = (0.36,0.36,0.36)
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').Transparency = 50


        self.form[5].progressBar.setValue(100) 


        # Switch to FEM GUI
        App.ActiveDocument.recompute()


        Gui.Control.closeDialog()
        Gui.activateWorkbench("FemWorkbench")
        FemGui.setActiveAnalysis(App.activeDocument().getObject(analysisName))

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)
        Gui.runCommand('Std_PerspectiveCamera',1)






        
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