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

import FreeCAD as App
import ObjectsFem



def genAnalysis(analysisName,constitutiveEQ):

    # See if analysis already exists
    try:
        test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(analysisName)[0] != None)
    except:
        test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(analysisName) != [])

    if test == False:
        # Analysis
        analysis_object = ObjectsFem.makeAnalysis(App.ActiveDocument,analysisName)

        # Solver 
        #solver_object = ObjectsFem.makeSolverCalculixCcxTools(App.ActiveDocument, "Project Chrono")
        #solver_object.GeometricalNonlinearity = 'linear'
        #solver_object.ThermoMechSteadyState = True
        #solver_object.MatrixSolverType = 'default'
        #solver_object.IterationsControlParameterTimeUse = False
        #analysis_object.addObject(solver_object)

        # Store Material
        material_object = ObjectsFem.makeMaterialSolid(App.ActiveDocument, constitutiveEQ)
        mat = material_object.Material
        mat['Name'] = constitutiveEQ
        material_object.Material = mat
        analysis_object.addObject(material_object)