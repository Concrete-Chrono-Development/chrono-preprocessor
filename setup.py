from setuptools import setup
import os
# from freecad.chronoWorkbench.version import __version__
# name: this is the name of the distribution.
# Packages using the same name here cannot be installed together

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
