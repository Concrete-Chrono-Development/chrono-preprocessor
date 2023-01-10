from setuptools import setup
import os
# from freecad.chronoConcrete.version import __version__
# name: this is the name of the distribution.
# Packages using the same name here cannot be installed together

version_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 
                            "freecad", "chronoConcrete", "version.py")
with open(version_path) as fp:
    exec(fp.read())

setup(name='freecad.chronoConcrete',
      version=str(__version__),
      packages=['freecad',
                'freecad.chronoConcrete'],
      maintainer="mtroemner",
      maintainer_email="mtroemner@gmail.com",
      url="TBD",
      description="template for a freecad extensions, installable with pip",
      install_requires=['numpy'], 
      include_package_data=True)
