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

from setuptools import setup
import os

version_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 
                            "freecad", "chronoWorkbench", "version.py")
with open(version_path) as fp:
    exec(fp.read())

setup(name='freecad.chronoWorkbench',
      version=str(__version__),
      packages=['freecad',
                'freecad.chronoWorkbench'],
      maintainer="mtroemner",
      maintainer_email="mtroemner@gmail.com",
      url="TBD",
      description="Chrono Workbench",
      install_requires=['numpy','math','time'], 
      include_package_data=True)
