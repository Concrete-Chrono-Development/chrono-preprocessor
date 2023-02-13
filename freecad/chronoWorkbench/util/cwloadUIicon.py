import os
from PySide import QtCore, QtGui
from freecad.chronoWorkbench import ICONPATH


def cwloadUIicon(form,icon):

    form.setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, icon))))