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

from freecad.chronoWorkbench                                     import ICONPATH
from freecad.chronoWorkbench                                     import GUIPATH
from freecad.chronoWorkbench                                     import TETGENPATH

from freecad.chronoWorkbench.util.cwloadUIfile                   import cwloadUIfile
from freecad.chronoWorkbench.util.cwloadUIicon                   import cwloadUIicon

from freecad.chronoWorkbench.generation.genAnalysis              import genAnalysis
from freecad.chronoWorkbench.generation.genGeometry              import genGeometry
from freecad.chronoWorkbench.generation.genSurfaceMesh           import genSurfMesh
from freecad.chronoWorkbench.generation.particleVol              import particleVol
from freecad.chronoWorkbench.generation.particleList             import particleList
from freecad.chronoWorkbench.generation.particleFaces            import particleFaces
from freecad.chronoWorkbench.generation.surfMeshSize             import surfMeshSize
from freecad.chronoWorkbench.generation.surfMeshExtents          import surfMeshExtents
from freecad.chronoWorkbench.generation.genParticle              import generateParticle
from freecad.chronoWorkbench.generation.genParticleMPI           import generateParticleMPI
from freecad.chronoWorkbench.generation.particleInsideCheck      import insideCheck
from freecad.chronoWorkbench.generation.particleOverlapCheck     import overlapCheck
from freecad.chronoWorkbench.generation.particleOverlapCheckMPI  import overlapCheckMPI
from freecad.chronoWorkbench.generation.readTetgen               import readTetgen
from freecad.chronoWorkbench.generation.genTetrahedralization    import genTetrahedralization
from freecad.chronoWorkbench.generation.genTesselation           import genTesselation
from freecad.chronoWorkbench.generation.genFacetData             import genFacetData

from freecad.chronoWorkbench.input.readInputsLDPM                import readInputs

from freecad.chronoWorkbench.output.mkVtkParticles               import mkVtkParticles
from freecad.chronoWorkbench.output.mkVtkFacets                  import mkVtkFacets
from freecad.chronoWorkbench.output.mkDataNodes                  import mkDataNodes
from freecad.chronoWorkbench.output.mkDataTets                   import mkDataTets
from freecad.chronoWorkbench.output.mkDataFacets                 import mkDataFacets



# Turn off error for divide by zero and invalid operations
np.seterr(divide='ignore', invalid='ignore')



#sys.executable = str(Path(App.ConfigGet('AppHomePath') + '/bin/python.exe'))
multiprocessing.set_executable(str(Path(App.ConfigGet('AppHomePath') + '/bin/pythonw.exe')))



class inputWindow_LDPM_CSL:
    def __init__(self):

        self.form = []

        # Load UI's for Side Panel
        self.form.append(cwloadUIfile("LDPM_CSL_modelProps.ui"))
        self.form.append(cwloadUIfile("LDPM_CSL_geometry.ui"))
        self.form.append(cwloadUIfile("LDPM_CSL_particles.ui"))        
        self.form.append(cwloadUIfile("LDPM_CSL_mixDesign.ui"))          
        self.form.append(cwloadUIfile("LDPM_CSL_additionalPara.ui"))       
        self.form.append(cwloadUIfile("LDPM_CSL_generation.ui"))

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
        QtCore.QObject.connect(self.form[5].generateLDPM, QtCore.SIGNAL("clicked()"), self.generateLDPM)
        QtCore.QObject.connect(self.form[5].generateCSL, QtCore.SIGNAL("clicked()"), self.generateCSL)
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
            outputDir] = readInputs(self.form)

        # Write parameters to file
        # Make output directory if does not exist
        outDir =  self.form[5].outputDir.text()
        try:
            os.mkdir(outDir)
        except:
            pass

        # Write parameters to file
        with open(Path(outDir + "/chronoWorkbench.cwPar"), "w") as f:
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
        self.form[0].matParaSet.setCurrentText(matParaSet)
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





    def generateLDPM(self):

        self.generation("LDPM")



    def generateCSL(self):

        self.generation("CSL")




    def generation(self, elementType):

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



        #with open(Path(tempPath + "file.txt"), "w") as f:
        #    path = "where python > " + str(Path(tempPath + "file.txt"))
        #    os.system(path)

        #with open(Path(tempPath + "file.txt"), "r") as f:
        #    output = f.readline()

        #new_string = output.replace("python", "pythonw")

        #print(new_string)

        



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
            outputDir] = readInputs(self.form)

        try:
            sieveCurveDiameter = ast.literal_eval(sieveCurveDiameter)
            sieveCurvePassing = ast.literal_eval(sieveCurvePassing)
        except:
            pass




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
        genGeo = genGeometry(dimensions,geoType,geoName,cadFile)
        self.form[5].progressBar.setValue(2) 

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)
        Gui.runCommand('Std_PerspectiveCamera',1)


        # Generate analysis objects
        self.form[5].statusWindow.setText("Status: Generating analysis objects.") 
        genAna = genAnalysis(analysisName,materialName)
        self.form[5].progressBar.setValue(3) 



        # Generate surface mesh
        self.form[5].statusWindow.setText("Status: Generating surface mesh.") 
        [vertices,edges,faces,tets] = genSurfMesh(analysisName,geoName,meshName,minPar,maxPar)
        self.form[5].progressBar.setValue(5) 



        # Gets extents of geometry
        [minC,maxC] = surfMeshExtents(vertices)

        self.form[5].statusWindow.setText("Status: Calculating input data.") 
        # Calculate required volume of particles and sieve curve data
        [parVolTotal,cdf,cdf1,kappa_i] = particleVol(wcRatio,airFrac,fullerCoef,cementC,cementDensity,densityWater,\
            flyashC,silicaC,scmC,flyashDensity,silicaDensity,scmDensity,fillerC,fillerDensity,\
            vertices,tets,minPar,maxPar,sieveCurveDiameter,sieveCurvePassing)

        # Temporary to skip over sieve curve option
        newSieveCurveD = 0
        NewSet = 0

        self.form[5].statusWindow.setText("Status: Calculating list of particles.") 
        # Calculate list of particle diameters for placement
        [maxParNum,parDiameterList] = particleList(parVolTotal,minPar,maxPar,newSieveCurveD,\
            cdf,kappa_i,NewSet,fullerCoef)

        # Calculation of surface mesh size
        maxEdgeLength = surfMeshSize(vertices,faces)

        # Generates points for all external triangles
        facePoints = particleFaces(vertices,faces)

        # Basic Calcs
        aggOffset = aggOffsetCoeff*minPar

        
        # Store coordinates of tets in new format
        coord1 = vertices[tets[:,0]-1]
        coord2 = vertices[tets[:,1]-1]
        coord3 = vertices[tets[:,2]-1]
        coord4 = vertices[tets[:,3]-1]




        verts = vertices[np.array(tets).flatten()-1]
        max_dist = np.max(np.sqrt(np.sum(verts**2, axis=1)))




        # Initialize empty particle nodes list outside geometry
        nodes = (np.zeros((len(parDiameterList),3))+2)*maxC




        self.form[5].statusWindow.setText('Status: Placing particles into geometry. (' + str(0) + '/' + str(len(parDiameterList)) + ')') 
        # Initialize values
        newMaxIter = 6
        particlesPlaced = 0

        




        if numCPU > 1:
        
            
            for increment in range(numIncrements-1):

                process_pool = multiprocessing.Pool(numCPU)

                outputMPI = process_pool.map(functools.partial(generateParticleMPI, facePoints,maxParNum, minC, maxC, vertices, \
                    tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                    maxPar,aggOffset,verbose,parDiameterList,maxEdgeLength,max_dist,nodes), parDiameterList[particlesPlaced:particlesPlaced+math.floor(len(parDiameterList)/numIncrements)])

                nodeMPI = np.array(outputMPI)[:,0:3]
                diameter = np.array(outputMPI)[:,3]
                newMaxIter = int(max(np.array(outputMPI)[:,4]))
                maxAttempts = int(max(np.array(outputMPI)[:,5]))

                particlesPlaced = particlesPlaced+len(np.array(outputMPI)[:,0:3])        

                for x in range(len(nodeMPI)):

                    # Store placed particles from this increment
                    nodes[particlesPlaced+x,:] = nodeMPI[x,:]

                    # Obtain extents for floating bin for node to test
                    binMin = np.array(([nodeMPI[x,0]-diameter[x]/2-maxPar/2-aggOffset,\
                        nodeMPI[x,1]-diameter[x]/2-maxPar/2-aggOffset,nodeMPI[x,2]-\
                        diameter[x]/2-maxPar/2-aggOffset]))
                    binMax = np.array(([nodeMPI[x,0]+diameter[x]/2+maxPar/2+aggOffset,\
                        nodeMPI[x,1]+diameter[x]/2+maxPar/2+aggOffset,nodeMPI[x,2]+\
                        diameter[x]/2+maxPar/2+aggOffset]))

                    # Check if particle overlapping any just added particles (ignore first one placed)
                    if x > 0:

                        overlap = overlapCheckMPI(nodeMPI[x,:],diameter[x],binMin,\
                            binMax,minPar,aggOffset,nodeMPI[0:x],diameter[0:x])

                        if overlap == True:

                            [newMaxIter,node,iterReq] = generateParticle(facePoints,\
                                parDiameterList[particlesPlaced+x], vertices, \
                                tets,newMaxIter,maxIter,minPar,\
                                maxPar,aggOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,nodes)
                            
                            nodes[particlesPlaced+x,:] = node[0,:]


                self.form[5].progressBar.setValue(95*((x)/len(parDiameterList))+6) 
                self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')


        # Generate particles for length of needed aggregate (not placed via MPI)
        for x in range(particlesPlaced,len(parDiameterList)):

            # Generate particle
            [newMaxIter,node,iterReq] = generateParticle(facePoints,parDiameterList[x],vertices,tets,newMaxIter,maxIter,minPar,maxPar,\
                aggOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,nodes)

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

            nodes[x,:] = node

        self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(len(parDiameterList)) + '/' + str(len(parDiameterList)) + ')')


        materialList = np.ones(len(parDiameterList))

        placementTime = round(time.time() - start_time,2)   
        nParticles = len(parDiameterList)

        # Create empty lists if not multi-material or cementStructure
        aggGrainsDiameterList, itzDiameterList, binderDiameterList, PoresDiameterList,\
            ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0,0,0,0









        tetTessTimeStart = time.time()


        self.form[5].statusWindow.setText("Status: Forming tetrahedralization.") 
        tetGen = genTetrahedralization(nodes,vertices,\
            faces,geoName,verbose,tempPath)
        self.form[5].progressBar.setValue(89) 



        [allNodes,allTets] = readTetgen(Path(tempPath + geoName \
        + '.node'),Path(tempPath + geoName + '.ele'))
        self.form[5].progressBar.setValue(90) 



        self.form[5].statusWindow.setText("Status: Forming tesselation.") 
        [tetFacets,facetCenters,facetAreas,facetNormals,tetn1,tetn2,tetPoints,allDiameters,facetPointData,facetCellData] = \
            genTesselation(allNodes,allTets,parDiameterList,minPar,\
            geoName)
        self.form[5].progressBar.setValue(95) 
        tetTessTime = round(time.time() - tetTessTimeStart,2)   




        # Store values for unused features
        edgeMaterialList = 0
        materialRule = 0
        multiMaterial = 'Off'
        cementStructure = 'Off'




        writeTimeStart = time.time()





     
        self.form[5].statusWindow.setText("Status: Generating facet data information.") 
        [facetData,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial] = genFacetData(\
            allNodes,allTets,tetFacets,facetCenters,facetAreas,facetNormals,tetn1,\
            tetn2,materialList,materialRule,multiMaterial,cementStructure,edgeMaterialList,facetCellData)
        self.form[5].progressBar.setValue(98) 







        self.form[5].statusWindow.setText("Status: Writing external facet data file.") 
        # Create file of external triangle facets for plotting of cells
        #externalFacetsFile = externalFacetFile(facetData,vertices,faces,geoName)





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


        print('Writing Node Data file.')

        mkDataNodes(geoName,tempPath,allNodes)

        print('Writing Tet Data file.')

        mkDataTets(geoName,tempPath,allTets)

        print('Writing Facet Data file.')

        # If data files requested, generate Facet File
        mkDataFacets(geoName,tempPath,facetData,facetPointData)



        print('Writing Particle Data file.')

        # If data files requested, generate Particle Data File
        #particleData(allNodes,allTets,parDiameterList,minPar,geoName)


        print('Writing visual files.')

        # If visuals requested, generate Particle VTK File
        mkVtkParticles(nodes,parDiameterList,materialList,geoName,tempPath)

        # If visuals requested, generate Facet VTK File
        mkVtkFacets(geoName,tempPath,facetPointData,facetCellData)



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
        Fem.insert(str(Path(outDir + outName + '/' + geoName + '-para-facet.000.vtk')),App.ActiveDocument.Name)        
        LDPMfacetsVTK = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facet_000')
        LDPMfacetsVTK.Label = 'LDPMfacetsVTK' 
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMfacetsVTK)
        LDPMfacetsVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-facet.000.vtk'))


        # Set visualization properties for facets
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facet_000').DisplayMode = u"Wireframe"
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facet_000').MaxFacesShowInner = 0
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facet_000').BackfaceCulling = False
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facet_000').ShapeColor = (0.36,0.36,0.36)


        # Set visualization properties for particle centers
        Gui.getDocument(App.ActiveDocument.Name).getObject(meshName).DisplayMode = u"Nodes"
        Gui.getDocument(App.ActiveDocument.Name).getObject(meshName).PointSize = 3.00
        Gui.getDocument(App.ActiveDocument.Name).getObject(meshName).PointColor = (0.00,0.00,0.00)










        # Store material properties

        # Property list from:
        #--- Smith, J., & Cusatis, G. (2017). Numerical analysis of projectile penetration and perforation of plain and fiber reinforced concrete slabs. 
        #--- International Journal for Numerical and Analytical Methods in Geomechanics, 41(3), 315-337.


        
        materialProps = [\
            "Density",\
            "CompressiveNormalModulus",\
            "PoissonsRatio",\
            "InitialHardeningModulusRatio",\
            "DensificationRatio",\
            "TransitionalStrainRatio",\
            "ShearStrengthRatio",\
            "InitialFriction",\
            "SofteningExponent",\
            "DeviatoricStrainThresholdRatio",\
            "DeviatoricDamageParameter",\
            "AsymptoticFriction",\
            "TransitionalStress",\
            "TensileStrength",\
            "TensileCharacteristicLength",\
            "ShearSofteningModulusRatio",\
            ]




        materialPropDesc = [\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            "Description coming soon...",\
            ]


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

        # Add appropriate material properties
        App.getDocument(App.ActiveDocument.Name).getObject(materialName).addProperty("App::PropertyString",'ConstitutiveEquationSet','Base','Set of constitutive equations.').ConstitutiveEquationSet=constitutiveEQ




        for x in range(len(materialProps)):
            App.getDocument(App.ActiveDocument.Name).getObject(materialName).addProperty("App::PropertyFloat",materialProps[x],elementType+" Parameters",materialPropDesc[x])#.Density=0.25


        # Add appropriate simulation properties
        for x in range(len(simProps)):
            App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyFloat",simProps[x],"Simulation",simPropDesc[x])#.Density=0.25
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyEnumeration","Solver","Simulation","Solver software").Solver=['Project Chrono']
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyEnumeration","IntegrationScheme","Simulation","Integrator type").IntegrationScheme=['Explicit']












        self.form[5].progressBar.setValue(100) 






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


class input_LDPM_CSL_Class():
    """My new command"""

    def GetResources(self):
        return {"Pixmap"  : os.path.join(ICONPATH, "ldpm.svg"), # the name of a svg file available in the resources
                "MenuText": "LDPM/CSL Generation",
                "ToolTip" : "Generation of an LDPM or CSL geometry"}

    def Activated(self):

        Gui.Control.showDialog(inputWindow_LDPM_CSL())

        return

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True

Gui.addCommand("mod_LDPM_CSL", input_LDPM_CSL_Class())