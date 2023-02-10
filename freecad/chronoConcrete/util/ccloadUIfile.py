import os
import FreeCADGui as Gui
from freecad.chronoConcrete import GUIPATH


def ccloadUIfile(filename):

        ui_file_A = os.path.join(GUIPATH,filename)
        a = Gui.PySideUic.loadUi(ui_file_A)

        return a