import os
from PySide import QtCore, QtGui
from freecad.chronoConcrete import ICONPATH


def ccloadUIicon(form,icon):

    form.setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, icon))))