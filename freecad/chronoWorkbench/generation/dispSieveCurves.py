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
## Developed by Northwestern University
## Primary Authors: Matthew Troemner
## ================================================================================
##
## This file contains the function to generate and display sieve curves for the 
## input sieve curve and generated particle size distribution.
##
## ================================================================================

# pyright: reportMissingImports=false
try:
    from FreeCAD.Plot import Plot
except ImportError:
    from freecad.plot import Plot
import numpy as np
import math
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as mticker

def dispSieveCurves(fullerCoef,sieveCurveDiameter,sieveCurvePassing,parDiameterList):

    """
    Variable List:
    --------------------------------------------------------------------------
    ### Inputs ###
    fullerCoef:              Fuller coefficient of the input particle size distribution
    sieveCurveDiameter:      List of diameters for the input sieve curve
    sieveCurvePassing:       List of percent passing for the input sieve curve
    parDiameterList:         List of diameters for the generated particle size distribution
    --------------------------------------------------------------------------
    ### Outputs ###
    None
    --------------------------------------------------------------------------
    """

    ############################################ TO DO ############################################

    # Calculations for sieve curve plotting for input sieve curve (discrete sieve case)
    
    # Calculations for sieve curve plotting for input sieve curve (fuller coefficient case)



    # Calculations for sieve curve plotting for generated particle size distribution
    totalVol = sum(4/3*math.pi*(parDiameterList/2)**3)
    passing = np.zeros(len(parDiameterList))
    for x in range(len(parDiameterList)):
        passing[x] = sum(4/3*math.pi*(parDiameterList[0:(len(parDiameterList)-x)]/2)**3)/totalVol*100
    

    # Generate plot of sieve curve
    Plot.figure("Particle Sieve Curve")
    Plot.plot(parDiameterList, passing, 'Simulated Data') 
 

    Plot.xlabel('Particle Diameter, $d$ (mm)') 
    Plot.ylabel('Percent Passing, $P$ (%)')

    Plot.grid(True)
    Plot.legend() 

