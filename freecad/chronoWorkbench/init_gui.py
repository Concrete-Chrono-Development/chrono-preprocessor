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

import os
import FreeCADGui as Gui #type: ignore
import FreeCAD as App #type: ignore
from PySide import QtGui #type: ignore
from FreeCADGui import Workbench #type: ignore

# Paths to Import
from freecad.chronoWorkbench import ICONPATH
from freecad.chronoWorkbench import GUIPATH

# Chrono Scripts to Import
from freecad.chronoWorkbench.modules import mod_LDPMCSL
from freecad.chronoWorkbench.modules import mod_LDPMCSL_gen
from freecad.chronoWorkbench.modules import mod_SPHDEM



class ChronoWorkbench(Gui.Workbench):
    """
    class which gets initiated at startup of the gui
    """

    MenuText = "Chrono Workbench"
    ToolTip = "A workbench for building LDPM, CSL, DEM, and SPH models for Project Chrono"
    Icon = os.path.join(ICONPATH, "ldpm.svg")
    toolbox = ["mod_LDPMCSL","mod_LDPMCSL_gen","mod_SPHDEM"] # a list of command names 


    def Initialize(self):
        """
        This function is called at the first activation of the workbench.
        here is the place to import all the commands
        """
        
        App.Console.PrintMessage("Switching to Chrono Workbench\n")
        App.Console.PrintMessage("A workbench for building LDPM, CSL, DEM, and SPH models for Project Chrono and other solvers.\n")

        self.appendToolbar("Tools", self.toolbox) # creates a new toolbar with your commands
        self.appendMenu("Tools", self.toolbox) # creates a new menu
        self.appendMenu(["Chrono Workbench"], self.toolbox) # appends a submenu to an existing menu



    def Activated(self):
        '''
        code which should be computed when a user switch to this workbench
        '''
        pass

    def Deactivated(self):
        '''
        code which should be computed when this workbench is deactivated
        '''
        pass

    def ContextMenu(self, recipient):
        """This function is executed whenever the user right-clicks on screen"""
        # "recipient" will be either "view" or "tree"
        self.appendContextMenu("My commands", self.list) # add commands to the context menu


    def GetClassName(self): 
        # This function is mandatory if this is a full Python workbench
        # This is not a template, the returned string should be exactly "Gui::PythonWorkbench"
        return "Gui::PythonWorkbench"



Gui.addWorkbench(ChronoWorkbench())



