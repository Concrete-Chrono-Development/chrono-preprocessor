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
import FreeCADGui as Gui
from freecad.chronoWorkbench import GUIPATH


def cwloadUIfile(filename):

        ui_file_A = os.path.join(GUIPATH,filename)
        a = Gui.PySideUic.loadUi(ui_file_A)

        return a