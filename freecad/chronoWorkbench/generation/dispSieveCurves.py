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

def dispSieveCurves(volFracPar, tetVolume, minPar, maxPar,fullerCoef,sieveCurveDiameter,sieveCurvePassing,parDiameterList):

    """
    Variable List:
    --------------------------------------------------------------------------
    ### Inputs ###
    volFracPar:              Volume fraction of particles in the geometry
    tetVolume:               Volume of the tetrahedral mesh
    minPar:                  Minimum particle diameter
    maxPar:                  Maximum particle diameter
    fullerCoef:              Fuller coefficient of the input particle size distribution
    sieveCurveDiameter:          List of diameters for the input sieve curve
    sieveCurvePassing:          List of percent passing for the input sieve curve
    parDiameterList:         List of diameters for the generated particle size distribution
    --------------------------------------------------------------------------
    ### Outputs ###
    Display of the input sieve curve and generated particle size distribution
    --------------------------------------------------------------------------
    """

    # Generate plot of sieve curve
    Plot.figure("Particle Sieve Curve")


    # Get volume of small particles and generated particles
    totalVol = sum(4/3*math.pi*(parDiameterList/2)**3)
    volParticles=volFracPar*tetVolume;
    volExtra=volParticles-totalVol;

    # Initialize Diameters
    parDiameterList = np.flip(parDiameterList)
    diameters = np.linspace(0,maxPar,num=1000)
    passingPercent=np.zeros(len(diameters))

    # Get Passing Percent of Placed Particles
    for x in range(len(diameters)):
        passing=parDiameterList[parDiameterList<diameters[x]]
        vol=sum(4/3*math.pi*(passing/2)**3)+volExtra
        passingPercent[x]=vol/volParticles*100









    # Calculations for sieve curve plotting for shifted generated particle size distribution (for comparison with Fuller Curve)
    if fullerCoef != 0:



        # Generate values for small particles
        passingPercent[0:500]=100*(diameters[0:500]/maxPar)**fullerCoef

        # Plotting
        Plot.plot(diameters[500:1000], passingPercent[500:1000], 'Simulated Data (Shifted)') 
        Plot.plot(diameters[0:500], passingPercent[0:500], 'Theoretical Curve') 

    else:


        Plot.plot(diameters[500:1000], passingPercent[500:1000], 'Simulated Data (Shifted)') 

        sieveCurveDiameter = np.asarray(sieveCurveDiameter, dtype=np.float32)
        sieveCurvePassing = np.asarray([x*100 for x in sieveCurvePassing], dtype=np.float32)
      
        Plot.plot(sieveCurveDiameter, sieveCurvePassing, 'Theoretical Curve (Adjusted)')

        ############################################ TO DO ############################################
        ####### Update to include plotting of discrete sieve curve    
    



    Plot.xlabel('Particle Diameter, $d$ (mm)') 
    Plot.ylabel('Percent Passing, $P$ (%)')

    Plot.grid(True)
    Plot.legend() 

