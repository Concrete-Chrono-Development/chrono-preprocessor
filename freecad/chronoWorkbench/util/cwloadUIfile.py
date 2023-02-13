import os
import FreeCADGui as Gui
from freecad.chronoWorkbench import GUIPATH


def cwloadUIfile(filename):

        ui_file_A = os.path.join(GUIPATH,filename)
        a = Gui.PySideUic.loadUi(ui_file_A)

        return a