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
## Function to assign icons from the 'icon' directory to modules in the left
## pane window of FreeCAD.
##
## ===========================================================================

import os
from PySide import QtCore, QtGui
from freecad.chronoWorkbench import ICONPATH


def cwloadUIicon(form,icon):

    form.setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, icon))))