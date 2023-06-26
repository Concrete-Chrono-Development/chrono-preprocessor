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

import os
import sys
import time
import tempfile
import numpy as np
from pathlib import Path
import multiprocessing
import functools
import math

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


from freecad.chronoWorkbench.util.cwloadUIfile                   import cwloadUIfile
from freecad.chronoWorkbench.util.cwloadUIicon                   import cwloadUIicon




# Turn off error for divide by zero and invalid operations
np.seterr(divide='ignore', invalid='ignore')



#sys.executable = str(Path(App.ConfigGet('AppHomePath') + '/bin/python.exe'))
multiprocessing.set_executable(str(Path(App.ConfigGet('AppHomePath') + '/bin/pythonw.exe')))



class inputWindow_SPHDEM:
    def __init__(self):

        self.form = []

        # Load UI's for Side Panel
        self.form.append(cwloadUIfile("ui_SPHDEM_geometry.ui"))
        self.form.append(cwloadUIfile("ui_SPHDEM_generation.ui"))

        # Label, Load Icons, and Initialize Panels
        #self.form[0].setWindowTitle("Model Settings")
        self.form[0].setWindowTitle("Geometry")
        #self.form[2].setWindowTitle("Particles")        
        #self.form[3].setWindowTitle("Mix Design")
        #self.form[4].setWindowTitle("Additional Parameters")
        self.form[1].setWindowTitle("Model Generation") 

        #cwloadUIicon(self.form[0],"FEM_MaterialMechanicalNonlinear.svg")
        cwloadUIicon(self.form[0],"PartDesign_AdditiveBox.svg")
        #cwloadUIicon(self.form[2],"Arch_Material_Group.svg")
        #cwloadUIicon(self.form[3],"FEM_ConstraintFlowVelocity.svg")
        #cwloadUIicon(self.form[4],"FEM_CreateNodesSet.svg")
        cwloadUIicon(self.form[1],"ldpm.svg")

        # Set initial output directory
        self.form[1].outputDir.setText(str(Path(App.ConfigGet('UserHomePath') + '/chronoWorkbench')))

        # Connect Open File Buttons
        #QtCore.QObject.connect(self.form[0].readFileButton, QtCore.SIGNAL("clicked()"), self.openFilePara)
        #QtCore.QObject.connect(self.form[0].readFileButton, QtCore.SIGNAL("clicked()"), self.openFileGeo)
        QtCore.QObject.connect(self.form[1].readDirButton, QtCore.SIGNAL("clicked()"), self.openDir)

        # Run generation for LDPM or CSL
        QtCore.QObject.connect(self.form[1].generateSPH, QtCore.SIGNAL("clicked()"), self.generation)
        QtCore.QObject.connect(self.form[1].generateDEM, QtCore.SIGNAL("clicked()"), self.generation)
        QtCore.QObject.connect(self.form[1].writePara, QtCore.SIGNAL("clicked()"), self.generation)



    def getStandardButtons(self):

        # Only show a close button
        # def accept() in no longer needed, since there is no OK button
        return int(QtGui.QDialogButtonBox.Close)

    def openFilePara(self):

        path = App.ConfigGet("UserHomePath")
        filetype = "CC Parameter input format (*.ccPar)"

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getOpenFileName(None,QString.fromLocal8Bit("Read a file parameter file"),path,             filetype)# type: ignore
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getOpenFileName(None, "Read a file parameter file", path,             filetype) #PySide
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        if OpenName == "":                                                            # if the name file are not selected then Abord process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            self.form[0].setupFile.setText(OpenName)

    def openFileGeo(self):

        path = App.ConfigGet("UserHomePath")
        filetype = "Supported formats (*.brep *.brp *.iges *.igs *.step *.stp);;\
                    BREP format       (*.brep *.brp);; \
                    IGES format       (*.iges *.igs);; \
                    STEP format       (*.step *.stp)"

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getOpenFileName(None,QString.fromLocal8Bit("Read a geometry file"),path,             filetype)# type: ignore
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





        elementType = "SPH"






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


        geoName = elementType + "geo" + str(0).zfill(3)
        meshName = elementType + "mesh" + str(0).zfill(3)
        analysisName = elementType + "analysis" + str(0).zfill(3)
        materialName = elementType + "material" + str(0).zfill(3)
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
            analysisName = elementType + "analysis" + str(i).zfill(3)
            materialName = elementType + "material" + str(i).zfill(3)
            dataFilesName = elementType + 'dataFiles'+ str(i).zfill(3)
            visualFilesName = elementType + 'visualFiles'+ str(i).zfill(3)
            try:
                test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0] != None)
            except:
                test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName) != [])

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


class input_SPHDEM_Class():
    """My new command"""

    def GetResources(self):
        return {"Pixmap"  : os.path.join(ICONPATH, "Arch_Material_Group.svg"), # the name of a svg file available in the resources
                "MenuText": "SPH/DEM Generation",
                "ToolTip" : "Generation of an SPH or DEM simulation"}

    def Activated(self):

        Gui.Control.showDialog(inputWindow_SPHDEM())

        return

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True

Gui.addCommand("mod_SPHDEM", input_SPHDEM_Class())