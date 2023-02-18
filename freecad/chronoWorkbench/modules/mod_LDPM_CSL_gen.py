## ================================================================================
## CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
##
## Copyright (c) 2023 
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be found
## in the LICENSE file at the top level of the distribution and at
## github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
##
## ================================================================================
## Author: Matthew Troemner
## ================================================================================
##
## Description coming soon...
##
##
## ================================================================================

import os
import sys
import time
import tempfile
import numpy as np
from pathlib import Path
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
import Spreadsheet

from PySide import QtCore, QtGui

from freecad.chronoWorkbench                                     import ICONPATH
from freecad.chronoWorkbench                                     import GUIPATH

from freecad.chronoWorkbench.util.cwloadUIfile                    import cwloadUIfile
from freecad.chronoWorkbench.util.cwloadUIicon                    import cwloadUIicon

from freecad.chronoWorkbench.output.mkChronoInput                  import mkChronoInput



class genWindow_LDPM_CSL:
    def __init__(self):

        self.form = []

        # Load UI's for Side Panel
        self.form.append(cwloadUIfile("LDPM_CSL_writeSimInfo.ui"))
        self.form.append(cwloadUIfile("LDPM_CSL_writeSimFiles.ui"))

        # Label, Load Icons, and Initialize Panels
        self.form[0].setWindowTitle("Confirm Output")
        self.form[1].setWindowTitle("Write File")

        cwloadUIicon(self.form[0],"ldpmOutput.svg")
        cwloadUIicon(self.form[1],"ldpmOutput.svg")


        # Set initial output directory
        self.form[1].outputDir.setText(str(Path(App.ConfigGet('UserHomePath') + '/chronoWorkbench')))

        # Connect Open File Buttons
        QtCore.QObject.connect(self.form[1].readDirButton, QtCore.SIGNAL("clicked()"), self.openDir)

        # Run generation for LDPM or CSL
        QtCore.QObject.connect(self.form[1].writeChrono, QtCore.SIGNAL("clicked()"), self.generation)



    def getStandardButtons(self):

        # Only show a close button
        # def accept() in no longer needed, since there is no OK button
        return int(QtGui.QDialogButtonBox.Close)



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



        print('Writing files.')

        elementType = "LDPM"
        materialProps = ["Alpha", "Normal Modulus", "etc."]
        materialPropsValues = (0,0,0)
        simProps = ["TBD", "TBD", "TBD", "TBD"]
        simPropsValues = (0,0,0,0)
        nodesFilename = "LDPMgeo000-data-nodes.dat"
        tetsFilename = "LDPMgeo000-data-tets.dat"
        facetsFilename = "LDPMgeo000-data-facets.dat"

        geoName = "TestGeo"
        geoType = "TestType"

        # Make output directory if does not exist
        outDir =  self.form[1].outputDir.text()
        try:
            os.mkdir(outDir)
        except:
            pass

        i = 0
        outName = '/' + geoName + geoType + str(i).zfill(3)
        while os.path.isdir(Path(outDir + outName)):
            i = i+1
            outName = '/' + geoName + geoType + str(i).zfill(3)

        try:
            os.mkdir(outDir + outName)
        except:
            pass


        # Make Project Chrono input file
        mkChronoInput(elementType, materialProps, materialPropsValues, simProps, simPropsValues, \
            nodesFilename, tetsFilename, facetsFilename, geoName, outDir, outName)







        
    # What to do when "Close" Button Clicked
    def reject(self):
         try:
             Gui.ActiveDocument.resetEdit()
             Gui.Control.closeDialog()
         except:
             Gui.Control.closeDialog()

        



class gen_LDPM_CSL_Class():
    """My new command"""

    def GetResources(self):
        return {"Pixmap"  : os.path.join(ICONPATH, "ldpmOutput.svg"), # the name of a svg file available in the resources
                "MenuText": "LDPM/CSL Simulation Inputs",
                "ToolTip" : "Generation of LDPM or CSL simulation inputs"}

    def Activated(self):

        Gui.Control.showDialog(genWindow_LDPM_CSL())

        return

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True

Gui.addCommand("mod_LDPM_CSL_gen", gen_LDPM_CSL_Class())