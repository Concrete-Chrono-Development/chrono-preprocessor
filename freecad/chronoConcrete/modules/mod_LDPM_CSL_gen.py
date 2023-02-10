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

from freecad.chronoConcrete                                     import ICONPATH
from freecad.chronoConcrete                                     import GUIPATH

from freecad.chronoConcrete.util.ccloadUIfile                    import ccloadUIfile
from freecad.chronoConcrete.util.ccloadUIicon                    import ccloadUIicon

from freecad.chronoConcrete.output.mkChronoInput                  import mkChronoInput



class genWindow_LDPM_CSL:
    def __init__(self):

        # Load UI's for Side Panel
        a = ccloadUIfile("LDPM_CSL_writeSimFiles.ui")
        self.form = [a]

        # Label, Load Icons, and Initialize Panels
        self.form[0].setWindowTitle("Write File")


        ccloadUIicon(self.form[0],"ldpmOutput.svg")


        # Set initial output directory
        self.form[0].outputDir.setText(str(Path(App.ConfigGet('UserHomePath') + '/chronoWorkbench')))

        # Connect Open File Buttons
        QtCore.QObject.connect(self.form[0].readDirButton, QtCore.SIGNAL("clicked()"), self.openDir)

        # Run generation for LDPM or CSL
        QtCore.QObject.connect(self.form[0].writeChrono, QtCore.SIGNAL("clicked()"), self.generation)



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
        materialProps = ("b" "h")
        materialPropsValues = (3,4)
        simProps = ("e" "s")
        simPropsValues = (1,2)
        nodesFilename = "test"
        tetsFilename = "test"
        
        geoName = "TestGeo"
        geoType = "TestType"

        # Make output directory if does not exist
        outDir =  self.form[0].outputDir.text()
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
            nodesFilename, tetsFilename, geoName, outDir, outName)







        
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