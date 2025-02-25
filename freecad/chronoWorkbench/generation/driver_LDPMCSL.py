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
## For U.S. Army ERDC Contract No. W9132T22C0015
## Primary Authors: Matthew Troemner
## ===========================================================================
##
## Description coming soon...
##
##
## ===========================================================================

# pyright: reportMissingImports=false

# Importing: standard
import os
import re
import shutil
import time
import tempfile
import numpy as np
from pathlib import Path
import multiprocessing
import functools
import math
import ast

# Importing: FreeCAD
import FreeCADGui as Gui
import FreeCAD as App
import Part
import Part,PartGui
import Mesh
import MeshPartGui, FreeCADGui
import MeshPart
import Mesh, Part, PartGui
import MaterialEditor
import ObjectsFem
import FemGui
import Fem
import femmesh.femmesh2mesh
from PySide import QtCore, QtGui


# Importing: generation
from freecad.chronoWorkbench.generation.calc_LDPMCSL_meshVolume           import calc_LDPMCSL_meshVolume
from freecad.chronoWorkbench.generation.calc_parVolume                    import calc_parVolume
from freecad.chronoWorkbench.generation.calc_sieveCurve                   import calc_sieveCurve
from freecad.chronoWorkbench.generation.calc_LDPMCSL_surfMeshSize         import calc_LDPMCSL_surfMeshSize
from freecad.chronoWorkbench.generation.calc_LDPMCSL_surfMeshExtents      import calc_LDPMCSL_surfMeshExtents
from freecad.chronoWorkbench.generation.check_particleOverlapMPI          import check_particleOverlapMPI
from freecad.chronoWorkbench.generation.check_multiMat_size               import check_multiMat_size
from freecad.chronoWorkbench.generation.check_multiMat_matVol             import check_multiMat_matVol
from freecad.chronoWorkbench.generation.gen_CSL_facetData                 import gen_CSL_facetData
from freecad.chronoWorkbench.generation.gen_LDPMCSL_tesselation           import gen_LDPMCSL_tesselation
from freecad.chronoWorkbench.generation.gen_LDPM_facetData                import gen_LDPM_facetData
from freecad.chronoWorkbench.generation.gen_LDPMCSL_analysis              import gen_LDPMCSL_analysis
from freecad.chronoWorkbench.generation.gen_LDPMCSL_facetfiberInt         import gen_LDPMCSL_facetfiberInt
from freecad.chronoWorkbench.generation.gen_LDPMCSL_fibers                import gen_LDPMCSL_fibers
from freecad.chronoWorkbench.generation.gen_LDPMCSL_flowEdges             import gen_LDPMCSL_flowEdges
from freecad.chronoWorkbench.generation.gen_LDPMCSL_geometry              import gen_LDPMCSL_geometry
from freecad.chronoWorkbench.generation.gen_LDPMCSL_initialMesh           import gen_LDPMCSL_initialMesh
from freecad.chronoWorkbench.generation.gen_particle                      import gen_particle
from freecad.chronoWorkbench.generation.gen_particleMPI                   import gen_particleMPI
from freecad.chronoWorkbench.generation.gen_particleList                  import gen_particleList
from freecad.chronoWorkbench.generation.gen_LDPMCSL_properties            import gen_LDPMCSL_properties
from freecad.chronoWorkbench.generation.gen_LDPMCSL_subParticle           import gen_LDPMCSL_subParticle
from freecad.chronoWorkbench.generation.gen_LDPMCSL_tetrahedralization    import gen_LDPMCSL_tetrahedralization
from freecad.chronoWorkbench.generation.gen_LDPMCSL_periodicMesh          import gen_LDPMCSL_periodicMesh
from freecad.chronoWorkbench.generation.gen_multiMat_refine               import gen_multiMat_refine
from freecad.chronoWorkbench.generation.gen_multiMat_reform               import gen_multiMat_reform
from freecad.chronoWorkbench.generation.gen_multiMat_assign               import gen_multiMat_assign
from freecad.chronoWorkbench.generation.sort_multiMat_voxels              import sort_multiMat_voxels
from freecad.chronoWorkbench.generation.sort_multiMat_mat                 import sort_multiMat_mat

# Importing: input
from freecad.chronoWorkbench.input.read_ctScan_file                       import read_ctScan_file
from freecad.chronoWorkbench.input.read_LDPMCSL_inputs                    import read_LDPMCSL_inputs
from freecad.chronoWorkbench.input.read_LDPMCSL_tetgen                    import read_LDPMCSL_tetgen
from freecad.chronoWorkbench.input.read_multiMat_file                     import read_multiMat_file

# Importing: output
from freecad.chronoWorkbench.output.mkVtk_particles                       import mkVtk_particles
from freecad.chronoWorkbench.output.mkVtk_LDPMCSL_facets                  import mkVtk_LDPMCSL_facets
from freecad.chronoWorkbench.output.mkVtk_LDPMCSL_fibers                  import mkVtk_LDPMCSL_fibers
from freecad.chronoWorkbench.output.mkVtk_LDPMCSL_flowEdges               import mkVtk_LDPMCSL_flowEdges
from freecad.chronoWorkbench.output.mkVtk_LDPMCSL_nonIntFibers            import mkVtk_LDPMCSL_nonIntFibers
from freecad.chronoWorkbench.output.mkVtk_LDPMCSL_projFacets              import mkVtk_LDPMCSL_projFacets
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleTetFacets            import mkVtk_LDPM_singleTetFacets
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleEdgeFacets           import mkVtk_LDPM_singleEdgeFacets
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleTetParticles         import mkVtk_LDPM_singleTetParticles
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleEdgeParticles        import mkVtk_LDPM_singleEdgeParticles
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleTet                  import mkVtk_LDPM_singleTet
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleEdge                 import mkVtk_LDPM_singleEdge
from freecad.chronoWorkbench.output.mkVtk_LDPM_singleCell                 import mkVtk_LDPM_singleCell
from freecad.chronoWorkbench.output.mkPy_LDPM_singleParaview              import mkPy_LDPM_singleParaview
from freecad.chronoWorkbench.output.mkPy_LDPM_singleParaviewLabels        import mkPy_LDPM_singleParaviewLabels
from freecad.chronoWorkbench.output.mkData_nodes                          import mkData_nodes
from freecad.chronoWorkbench.output.mkData_LDPMCSL_tets                   import mkData_LDPMCSL_tets
from freecad.chronoWorkbench.output.mkData_LDPMCSL_edges                  import mkData_LDPMCSL_edges
from freecad.chronoWorkbench.output.mkData_LDPMCSL_facets                 import mkData_LDPMCSL_facets
from freecad.chronoWorkbench.output.mkData_LDPMCSL_facetfiberInt          import mkData_LDPMCSL_facetfiberInt
from freecad.chronoWorkbench.output.mkData_LDPMCSL_facetsVertices         import mkData_LDPMCSL_facetsVertices
from freecad.chronoWorkbench.output.mkData_LDPMCSL_faceFacets             import mkData_LDPMCSL_faceFacets
from freecad.chronoWorkbench.output.mkData_LDPMCSL_flowEdges              import mkData_LDPMCSL_flowEdges
from freecad.chronoWorkbench.output.mkData_particles                      import mkData_particles
from freecad.chronoWorkbench.output.mkDisp_sieveCurves                    import mkDisp_sieveCurves
from freecad.chronoWorkbench.output.mkIges_LDPMCSL_flowEdges              import mkIges_LDPMCSL_flowEdges



def driver_LDPMCSL(self,fastGen,tempPath):

    # Read in inputs from input panel
    [setupFile, constitutiveEQ, matParaSet, \
        numCPU, numIncrements,maxIter,placementAlg,\
        geoType, dimensions, cadFile,\
        minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
        wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
        cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, \
        fillerC, fillerDensity, airFrac2,\
        htcToggle, htcLength,\
        fiberToggle, fiberCutting, fiberDiameter, fiberLength, fiberVol, fiberOrientation1, fiberOrientation2, fiberOrientation3, fiberPref, fiberFile, fiberIntersections,\
        multiMatToggle,aggFile,multiMatFile,multiMatRule,\
        grainAggMin, grainAggMax, grainAggFuller, grainAggSieveD, grainAggSieveP,\
        grainITZMin, grainITZMax, grainITZFuller, grainITZSieveD, grainITZSieveP,\
        grainBinderMin, grainBinderMax, grainBinderFuller, grainBinderSieveD, grainBinderSieveP,\
        periodicToggle,\
        outDir, dataFilesGen, visFilesGen, singleTetGen, modelType] = read_LDPMCSL_inputs(self.form)

    # Make output directory if does not exist
    try:
        os.mkdir(outDir)
    except:
        pass

    i = 0
    # Use single names for geoTypes
    if geoType in ["Box","Cylinder","Cone","Sphere","Ellipsoid","Prism","Dogbone","Custom"]:
        geoTypeOutName = geoType
    elif geoType == "Notched Prism - Semi Circle":
        geoTypeOutName = "NotchedPrismSemiCircle"
    elif geoType == "Notched Prism - Square":
        geoTypeOutName = "NotchedPrismSquare"
    elif geoType == "Notched Prism - Ellipse":
        geoTypeOutName = "NotchedPrismEllipse"
    elif geoType == "Import CAD or Mesh":
        geoTypeOutName = "ImportedFile"

    if modelType in ["Confinement Shear Lattice (CSL) - LDPM Style ",\
                        "Confinement Shear Lattice (CSL) - Original"]:
        elementType = "CSL"
    else:
        elementType = "LDPM"

    geoName = elementType + "geo" + str(0).zfill(3)

    outName = '/' + geoName + geoTypeOutName + str(i).zfill(3)
    while os.path.isdir(Path(outDir + outName)):
        i = i+1
        outName = '/' + geoName + geoTypeOutName + str(i).zfill(3)

    # Move existing files to selected output directory and remake temp directory
    print('Moving files.')    
    shutil.move(tempPath, outDir + outName)
    try:
        os.mkdir(tempPath)
    except:
        pass


    # Initialize code start time to measure performance
    start_time = time.time()


    # Store document
    docGui = Gui.activeDocument()

    # Make new document and set view if does not exisit
    try:
        docGui.activeView().viewAxonometric()
    except:
        App.newDocument("Unnamed")
        docGui = Gui.activeDocument()
        docGui.activeView().viewAxonometric()
    Gui.runCommand('Std_PerspectiveCamera',1)

    try:
        sieveCurveDiameter = ast.literal_eval(sieveCurveDiameter)
        sieveCurvePassing = ast.literal_eval(sieveCurvePassing)
    except:
        pass

    try:
        grainAggSieveD = ast.literal_eval(grainAggSieveD)
        grainAggSieveP = ast.literal_eval(grainAggSieveP)
    except:
        pass

    try:
        grainITZSieveD = ast.literal_eval(grainITZSieveD)
        grainITZSieveP = ast.literal_eval(grainITZSieveP)
    except:
        pass


    try:
        grainBinderSieveD = ast.literal_eval(grainBinderSieveD)
        grainBinderSieveP = ast.literal_eval(grainBinderSieveP)
    except:
        pass






    if fillerC > 0:
        airFrac = airFrac2
    else:
        airFrac = airFrac1
    
    parOffsetCoeff = 0.2                                    # Minimum distance between particles factor 
    verbose = "On"

    self.form[5].progressBar.setValue(1) 
    self.form[5].statusWindow.setText("Status: Generating objects.") 




    meshName = elementType + "mesh" + str(0).zfill(3)
    analysisName = elementType + "analysis"
    materialName = elementType + "material"
    dataFilesName = elementType + 'dataFiles'+ str(0).zfill(3)
    visualFilesName = elementType + 'visualFiles'+ str(0).zfill(3)

    i = 0
    try:
        test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0] != None)
    except:
        test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName) != [])

    while test == True:
        i = i+1
        geoName = elementType + "geo" + str(i).zfill(3)
        meshName = elementType + "mesh" + str(i).zfill(3)
        dataFilesName = elementType + 'dataFiles'+ str(i).zfill(3)
        visualFilesName = elementType + 'visualFiles'+ str(i).zfill(3)
        try:
            test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0] != None)
        except:
            test = (App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName) != [])

    # Generate geometry
    self.form[5].statusWindow.setText("Status: Generating geometry.") 
    genGeo = gen_LDPMCSL_geometry(dimensions,geoType,geoName,cadFile)
    self.form[5].progressBar.setValue(2) 

    # Set view
    docGui.activeView().viewAxonometric()
    Gui.SendMsgToActiveView("ViewFit")
    Gui.runCommand('Std_DrawStyle',6)
    Gui.runCommand('Std_PerspectiveCamera',1)


    # Generate analysis objects
    self.form[5].statusWindow.setText("Status: Generating analysis objects.") 
    genAna = gen_LDPMCSL_analysis(analysisName,materialName)
    self.form[5].progressBar.setValue(3) 


    # Generate surface mesh
    self.form[5].statusWindow.setText("Status: Generating surface mesh.") 
    if periodicToggle == "On":
        geoType = 'Import CAD or Mesh' # This is a hack to improve the visualization, as the periodic meshing is a quasi-imported mesh
        cadFile = gen_LDPMCSL_periodicMesh(cadFile,analysisName,geoName,meshName,minPar,dimensions,tempPath)
        [meshVertices,meshTets,surfaceNodes,surfaceFaces] = gen_LDPMCSL_initialMesh(cadFile,analysisName,geoName,meshName,minPar)
    else:
        [meshVertices,meshTets,surfaceNodes,surfaceFaces] = gen_LDPMCSL_initialMesh(cadFile,analysisName,geoName,meshName,minPar)

    self.form[5].progressBar.setValue(5) 






    # Gets extents of geometry
    [minC,maxC] = calc_LDPMCSL_surfMeshExtents(meshVertices)







    # Convert density to Kg/m3
    cementC = cementC * (1.0E+12)
    flyashC = flyashC * (1.0E+12)
    silicaC = silicaC * (1.0E+12)
    scmC = scmC * (1.0E+12)
    fillerC = fillerC * (1.0E+12)
    cementDensity = cementDensity * (1.0E+12)
    flyashDensity = flyashDensity * (1.0E+12)
    silicaDensity = silicaDensity * (1.0E+12)
    scmDensity = scmDensity * (1.0E+12)
    fillerDensity = fillerDensity * (1.0E+12)
    densityWater = densityWater * (1.0E+12)




    self.form[5].statusWindow.setText("Status: Calculating input data.") 
    

    # Gets volume of geometry
    tetVolume = calc_LDPMCSL_meshVolume(meshVertices,meshTets)





    # Calculation of surface mesh size
    maxEdgeLength = calc_LDPMCSL_surfMeshSize(meshVertices,surfaceFaces)


    # Basic Calcs
    parOffset = parOffsetCoeff*minPar

    
    # Store coordinates of meshTets in new format
    coord1 = meshVertices[meshTets[:,0]-1]
    coord2 = meshVertices[meshTets[:,1]-1]
    coord3 = meshVertices[meshTets[:,2]-1]
    coord4 = meshVertices[meshTets[:,3]-1]



    verts = meshVertices[np.array(meshTets).flatten()-1]
    max_dist = np.max(np.sqrt(np.sum(verts**2, axis=1)))




    if fastGen == True:

        ########################## Alternative Route to Farm Out Particle Processes ##############################

            # Make a temporary file that will be used to store the parameters and then run the generation:

    # Write these seven matrices to temporary files that will later be read back in by the generation script:
    # coord1, coord2, coord3, coord4, meshVertices, meshTets, surfaceNodes
        
        np.save(tempPath + "coord1.npy", coord1)
        np.save(tempPath + "coord2.npy", coord2)
        np.save(tempPath + "coord3.npy", coord3)
        np.save(tempPath + "coord4.npy", coord4)
        np.save(tempPath + "meshVertices.npy", meshVertices)
        np.save(tempPath + "meshTets.npy", meshTets)
        np.save(tempPath + "surfaceNodes.npy", surfaceNodes)

        # Get the current directory 
        currentDir = os.path.dirname(os.path.realpath(__file__))


        with open(Path(currentDir + "/tempGen.py"), "w") as f:
            f.write("""\n
# ================================================================================
# CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
#     
# ================================================================================
# Chrono Workbench Parameter File
# ================================================================================
#
# Chrono Workbench developed by Northwestern University
#
# ================================================================================

from gen_LDPMCSL_multiStep   import gen_LDPMCSL_multiStep                     
                    
            \n\n""")
            f.write('tempPath = r"' + tempPath + '"\n')
            f.write('numCPU = ' + str(numCPU) + "\n")
            f.write('numIncrements = ' + str(numIncrements) + "\n")
            f.write('maxIter = ' + str(maxIter) + "\n")
            f.write('parOffset = ' + str(parOffset) + "\n")
            f.write('maxEdgeLength = ' + str(maxEdgeLength) + "\n")
            f.write('max_dist = ' + str(max_dist) + "\n")
            f.write('minPar = ' + str(minPar) + "\n")
            f.write('maxPar = ' + str(maxPar) + "\n")
            if fullerCoef == "":
                f.write("fullerCoef = None\n")
            else:
                f.write("fullerCoef = " + str(fullerCoef) + "\n")
            if sieveCurveDiameter == "":
                f.write('sieveCurveDiameter = ""\n')
            else:
                f.write("sieveCurveDiameter = " + str(sieveCurveDiameter) + "\n")
            if sieveCurvePassing == "":
                f.write("sieveCurvePassing = None\n")
            else:
                f.write("sieveCurvePassing = " + str(sieveCurvePassing) + "\n")
            f.write("wcRatio = " + str(wcRatio) + "\n")
            f.write("densityWater = " + str(densityWater) + "\n")
            f.write("cementC = " + str(cementC) + "\n")
            f.write("flyashC = " + str(flyashC) + "\n")
            f.write("silicaC = " + str(silicaC) + "\n")
            f.write("scmC = " + str(scmC) + "\n")
            f.write("cementDensity = " + str(cementDensity) + "\n")
            f.write("flyashDensity = " + str(flyashDensity) + "\n")
            f.write("silicaDensity = " + str(silicaDensity) + "\n")
            f.write("scmDensity = " + str(scmDensity) + "\n")
            f.write("airFrac = " + str(airFrac) + "\n")
            f.write("fillerC = " + str(fillerC) + "\n")
            f.write("fillerDensity = " + str(fillerDensity) + "\n")
            f.write("tetVolume = " + str(tetVolume) + "\n")
            f.write("minC = [" + str(minC[0]) + ", " + str(minC[1]) + ", " + str(minC[2]) + "]\n")
            f.write("maxC = [" + str(maxC[0]) + ", " + str(maxC[1]) + ", " + str(maxC[2]) + "]\n")
            f.write('verbose = "' + str(verbose) + '"\n')
            if multiMatToggle == "On":
                f.write('multiMatToggle = "' + multiMatToggle + '"\n')
                f.write('multiMatFile = "' + multiMatFile + '"\n')
                f.write('aggFile = "' + aggFile + '"\n')
                f.write('multiMatRule = ' + str(multiMatRule) + '\n')
                f.write("grainAggMin = " + str(grainAggMin) + "\n")
                f.write("grainAggMax = " + str(grainAggMax) + "\n")
                if grainAggFuller == "":
                    f.write("grainAggFuller = None\n")
                else:
                    f.write("grainAggFuller = " + str(grainAggFuller) + "\n")
                if grainAggSieveD == "":
                    f.write('grainAggSieveD = ""\n')
                else:
                    f.write("grainAggSieveD = " + str(grainAggSieveD) + "\n")
                if grainAggSieveP == "":
                    f.write("grainAggSieveP = None\n")
                else:
                    f.write("grainAggSieveP = " + str(grainAggSieveP) + "\n")
                f.write("grainITZMin = " + str(grainITZMin) + "\n")
                f.write("grainITZMax = " + str(grainITZMax) + "\n")
                if grainITZFuller == "":
                    f.write("grainITZFuller = None\n")
                else:
                    f.write("grainITZFuller = " + str(grainITZFuller) + "\n")
                if grainITZSieveD == "":
                    f.write('grainITZSieveD = ""\n')  
                else:
                    f.write("grainITZSieveD = " + str(grainITZSieveD) + "\n")
                if grainITZSieveP == "":
                    f.write("grainITZSieveP = None\n")
                else:
                    f.write("grainITZSieveP = " + str(grainITZSieveP) + "\n")
                f.write("grainBinderMin = " + str(grainBinderMin) + "\n")
                f.write("grainBinderMax = " + str(grainBinderMax) + "\n")
                if grainBinderFuller == "":
                    f.write("grainBinderFuller = None\n")
                else:
                    f.write("grainBinderFuller = " + str(grainBinderFuller) + "\n")
                if grainBinderSieveD == "":
                    f.write('grainBinderSieveD = ""\n')
                else:
                    f.write("grainBinderSieveD = " + str(grainBinderSieveD) + "\n")
                if grainBinderSieveP == "":
                    f.write("grainBinderSieveP = None\n")
                else:
                    f.write("grainBinderSieveP = " + str(grainBinderSieveP) + "\n")
            else:
                f.write('multiMatToggle = "' + multiMatToggle + '"\n')
                f.write('aggFile = None\n')
                f.write('multiMatFile = None\n')
                f.write('multiMatRule = None\n')
                f.write("grainAggMin = None\n")
                f.write("grainAggMax = None\n")
                f.write("grainAggFuller = None\n")
                f.write("grainAggSieveD = None\n")
                f.write("grainAggSieveP = None\n")
                f.write("grainITZMin = None\n")
                f.write("grainITZMax = None\n")
                f.write("grainITZFuller = None\n")
                f.write("grainITZSieveD = None\n")
                f.write("grainITZSieveP = None\n")
                f.write("grainBinderMin = None\n")
                f.write("grainBinderMax = None\n")
                f.write("grainBinderFuller = None\n")
                f.write("grainBinderSieveD = None\n")
                f.write("grainBinderSieveP = None\n")


            f.write("""

def main():
                
    generation = gen_LDPMCSL_multiStep(tempPath, numCPU, numIncrements, maxIter, parOffset, maxEdgeLength, max_dist, minPar, maxPar, sieveCurveDiameter, sieveCurvePassing, wcRatio, cementC, airFrac, fullerCoef, flyashC, silicaC, scmC, fillerC, flyashDensity, silicaDensity, scmDensity, fillerDensity, cementDensity, densityWater, multiMatToggle, aggFile, multiMatFile, grainAggMin, grainAggMax, grainAggFuller, grainAggSieveD, grainAggSieveP, grainBinderMin, grainBinderMax, grainBinderFuller, grainBinderSieveD, grainBinderSieveP, grainITZMin, grainITZMax, grainITZFuller, grainITZSieveD, grainITZSieveP, tetVolume, minC, maxC, verbose)
                
                
if __name__ == '__main__':
    main()
            """)

        
        
        
        # Run the generation   
        os.system("python " + str(Path(currentDir + "/tempGen.py")))

        # Read the temporary internalNodes file
        internalNodes = np.load(tempPath + "internalNodes.npy")

        # Read the temporary materialList file
        materialList = np.load(tempPath + "materialList.npy")

        # Read the temporary parDiameterList file
        parDiameterList = np.load(tempPath + "parDiameterList.npy")

        # Read the particleID file
        particleID = np.load(tempPath + "particleID.npy")

        if multiMatToggle == "Off":
            # Read the volFracPar file
            volFracPar = np.load(tempPath + "volFracPar.npy")

        # Remove the temporary files
        os.remove(Path(currentDir + "/tempGen.py"))

        # Remove the temporary data files
        os.remove(tempPath + "internalNodes.npy")
        os.remove(tempPath + "materialList.npy")
        os.remove(tempPath + "parDiameterList.npy")
        if multiMatToggle == "Off":
            os.remove(tempPath + "volFracPar.npy")
        os.remove(tempPath + "coord1.npy")
        os.remove(tempPath + "coord2.npy")
        os.remove(tempPath + "coord3.npy")
        os.remove(tempPath + "coord4.npy")
        os.remove(tempPath + "meshVertices.npy")
        os.remove(tempPath + "meshTets.npy")
        os.remove(tempPath + "surfaceNodes.npy")
        os.remove(tempPath + "particleID.npy")



    else:





        #################### Begin Setting Up Particles and Materials (Normal Method) ##############################

        if multiMatToggle == "On":


            # Read in aggregate file
            try:
                [multiMatX,multiMatY,multiMatZ,multiMatRes,aggDistinctVoxels] = read_multiMat_file(aggFile)
            except:
                pass


            # Read in multi-material file
            [multiMatX,multiMatY,multiMatZ,multiMatRes,multiMatVoxels] = read_multiMat_file(multiMatFile)


            # Confirm if the voxelated multi-material file is larger than the provided geometry
            topoCheck = check_multiMat_size(multiMatX,multiMatY,multiMatZ,multiMatRes,minC,maxC)


            # Organize and store voxels of each material
            [aggVoxels,itzVoxels,binderVoxels,aggVoxelIDs] = sort_multiMat_voxels(multiMatVoxels)

            # Organize and store voxels of aggregate with the distinct ID file
            try:
                [aggVoxels,discard2,discard3,aggVoxelIDs] = sort_multiMat_voxels(aggDistinctVoxels)
            except:
                pass



            # Do calculations for aggregate, binder, and ITZ
            for i in range(3):
                
                if i == 0:
                    [grainMin,grainMax,grainFuller,grainSieveD,grainSieveP] = [grainAggMin,grainAggMax,grainAggFuller,grainAggSieveD,grainAggSieveP]
                elif i == 1:
                    [grainMin,grainMax,grainFuller,grainSieveD,grainSieveP] = [grainBinderMin,grainBinderMax,grainBinderFuller,grainBinderSieveD,grainBinderSieveP]
                elif i == 2:
                    [grainMin,grainMax,grainFuller,grainSieveD,grainSieveP] = [grainITZMin,grainITZMax,grainITZFuller,grainITZSieveD,grainITZSieveP]


                # Shift sieve curve if needed
                if grainSieveD != (0 or None or [] or ""):
                    [newGrainSieveCurveD,newGrainSieveCurveP,grainNewSet,grainW_min,grainW_max] = calc_sieveCurve(grainMin,grainMax,grainSieveD,grainSieveP)
                else:
                    newGrainSieveCurveD,newGrainSieveCurveP,grainNewSet,grainW_min,grainW_max = 0, 0, 0, 0, 0

                # Calculates volume of each set of grains
                [volGrainFracPar,volGrains,cdf,cdf1,kappa_i] = calc_parVolume(tetVolume*len(aggVoxels)/(len(aggVoxels)+len(itzVoxels)+len(binderVoxels)), wcRatio, cementC,
                                                            airFrac, grainFuller, 
                                                            flyashC, silicaC, scmC, fillerC,
                                                            flyashDensity, silicaDensity, 
                                                            scmDensity, fillerDensity, cementDensity,
                                                            densityWater, grainMin, grainMax,
                                                            newGrainSieveCurveD, newGrainSieveCurveP, 
                                                            grainNewSet, grainW_min, grainW_max)

                # Generates list of needed grains
                [maxGrainsNum,grainsDiameterList] = gen_particleList(volGrains,grainMin,grainMax,newGrainSieveCurveD,cdf,kappa_i,grainNewSet,grainFuller)

                if i == 0:
                    aggGrainsDiameterList = grainsDiameterList
                elif i == 1:
                    binderGrainsDiameterList = grainsDiameterList
                elif i == 2:
                    itzGrainsDiameterList = grainsDiameterList

            # Combine all grain lists (in order of aggregate > ITZ > binder -- the order they will be placed in the geometry)
            parDiameterList = np.concatenate((aggGrainsDiameterList,itzGrainsDiameterList,binderGrainsDiameterList))

            # Initialize empty list of all nodes outside geometry
            internalNodes = (np.zeros((len(aggGrainsDiameterList)+\
                len(binderGrainsDiameterList)+len(itzGrainsDiameterList),3))+2)*maxC




        if multiMatToggle == "Off":


            # Shift sieve curve if needed
            if sieveCurveDiameter != (0 or None or [] or ""):
                # Shifts sieve curve to appropriate range
                [newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max] = calc_sieveCurve(minPar, maxPar, sieveCurveDiameter, sieveCurvePassing)
            else:
                newSieveCurveD, newSieveCurveP, w_min, w_max, NewSet = 0, 0, 0, 0, 0

            # Calculates volume of particles needed
            [volFracPar, parVolTotal, cdf, cdf1, kappa_i] = calc_parVolume(tetVolume, wcRatio, cementC,
                                                        airFrac, fullerCoef, 
                                                        flyashC, silicaC, scmC, fillerC,
                                                        flyashDensity, silicaDensity, 
                                                        scmDensity, fillerDensity, cementDensity,
                                                        densityWater, minPar, maxPar,
                                                        newSieveCurveD, newSieveCurveP, 
                                                        NewSet, w_min, w_max)



            self.form[5].statusWindow.setText("Status: Calculating list of particles.") 
            # Calculate list of particle diameters for placement
            [maxParNum,parDiameterList] = gen_particleList(parVolTotal,minPar,maxPar,newSieveCurveD,\
                cdf,kappa_i,NewSet,fullerCoef)
        
            # Initialize empty particle nodes list outside geometry
            internalNodes = (np.zeros((len(parDiameterList),3))+2)*maxC
        












        ########################## Begin Placing Particles ##############################

        


        self.form[5].statusWindow.setText('Status: Placing particles into geometry. (' + str(0) + '/' + str(len(internalNodes)) + ')') 
        
        # Initialize values
        newMaxIter = 6
        particlesPlaced = 0


        # Initialize particleID list of length of internalNodes
        particleID = np.zeros(len(internalNodes))
        
        if multiMatToggle == "On":

            if len(itzVoxels) > 0:
                for i in range(3):
                    # Place in order of aggregate > ITZ > binder
                    if i == 0:
                        [grainsDiameterList,voxels,grainMin,grainMax,voxelIDs] = [aggGrainsDiameterList,aggVoxels,grainAggMin,grainAggMax,aggVoxelIDs]
                    elif i == 1:
                        [grainsDiameterList,voxels,grainMin,grainMax,voxelIDs] = [itzGrainsDiameterList,itzVoxels,grainITZMin,grainITZMax,0]
                    elif i == 2:
                        [grainsDiameterList,voxels,grainMin,grainMax,voxelIDs] = [binderGrainsDiameterList,binderVoxels,grainBinderMin,grainBinderMax,0]

                    # Generate particles for length of needed aggregate (not placed via MPI)
                    for x in range(particlesPlaced,len(grainsDiameterList)):

                        # Generate particle
                        [newMaxIter,node,iterReq,particleID[x]] = gen_LDPMCSL_subParticle(surfaceNodes,grainsDiameterList[x],meshVertices,meshTets,newMaxIter,maxIter,grainMin,grainMax,\
                            parOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes,\
                            multiMatX,multiMatY,multiMatZ,multiMatRes,voxels,voxelIDs,minC,maxC)

                        # Update progress bar every 1% of placement
                        if x % np.rint(len(grainsDiameterList)/100) == 0:
                            self.form[5].progressBar.setValue(80*((x)/len(grainsDiameterList))+6) 

                        if len(grainsDiameterList)<=1000:
                            # Update number particles placed every 1%
                            if x % np.rint(len(grainsDiameterList)/100) == 0:
                                self.form[5].statusWindow.setText("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')
                        elif len(grainsDiameterList)<=10000:
                            # Update number particles placed every 0.1%
                            if x % np.rint(len(grainsDiameterList)/1000) == 0:
                                self.form[5].statusWindow.setText("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')
                        else:
                            # Update number particles placed every 0.01%
                            if x % np.rint(len(grainsDiameterList)/10000) == 0:
                                self.form[5].statusWindow.setText("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')

                        if i == 0:
                            internalNodes[x,:] = node
                        elif i == 1:
                            internalNodes[x+len(aggGrainsDiameterList),:] = node
                        elif i == 2:
                            internalNodes[x+len(aggGrainsDiameterList)+len(itzGrainsDiameterList),:] = node
            


                    self.form[5].statusWindow.setText("Status: Placing material " + str(i) + " grains into geometry. (" + str(len(grainsDiameterList)) + '/' + str(len(grainsDiameterList)) + ')')

                materialList = np.concatenate((np.ones(len(aggGrainsDiameterList))*3,np.ones(len(itzGrainsDiameterList))*1, np.ones(len(binderGrainsDiameterList))*2))
            
                print(materialList)

                # Set minimum particle to be smallest of the three materials 
                minPar = min(grainAggMin,grainITZMin,grainBinderMin)



            else: 
                for i in range(2):
                    # Place in order of aggregate > binder
                    if i == 0:
                        [grainsDiameterList,voxels,grainMin,grainMax,voxelIDs] = [aggGrainsDiameterList,aggVoxels,grainAggMin,grainAggMax,aggVoxelIDs]
                    elif i == 1:
                        [grainsDiameterList,voxels,grainMin,grainMax,voxelIDs] = [binderGrainsDiameterList,binderVoxels,grainBinderMin,grainBinderMax,0]

                    # Generate particles for length of needed aggregate (not placed via MPI)
                    for x in range(particlesPlaced,len(grainsDiameterList)):

                        # Generate particle
                        [newMaxIter,node,iterReq,particleID[x]] = gen_LDPMCSL_subParticle(surfaceNodes,grainsDiameterList[x],meshVertices,meshTets,newMaxIter,maxIter,grainMin,grainMax,\
                            parOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes,\
                            multiMatX,multiMatY,multiMatZ,multiMatRes,voxels,voxelIDs,minC,maxC)

                        # Update progress bar every 1% of placement
                        if x % np.rint(len(grainsDiameterList)/100) == 0:
                            self.form[5].progressBar.setValue(80*((x)/len(grainsDiameterList))+6) 

                        if len(grainsDiameterList)<=1000:
                            # Update number particles placed every 1%
                            if x % np.rint(len(grainsDiameterList)/100) == 0:
                                self.form[5].statusWindow.setText("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')
                        elif len(grainsDiameterList)<=10000:
                            # Update number particles placed every 0.1%
                            if x % np.rint(len(grainsDiameterList)/1000) == 0:
                                self.form[5].statusWindow.setText("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')
                        else:
                            # Update number particles placed every 0.01%
                            if x % np.rint(len(grainsDiameterList)/10000) == 0:
                                self.form[5].statusWindow.setText("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')

                        if i == 0:
                            internalNodes[x,:] = node
                        elif i == 1:
                            internalNodes[x+len(aggGrainsDiameterList),:] = node
            

                    self.form[5].statusWindow.setText("Status: Placing material " + str(i) + " grains into geometry. (" + str(len(grainsDiameterList)) + '/' + str(len(grainsDiameterList)) + ')')

                materialList = np.concatenate((np.ones(len(aggGrainsDiameterList))*3,np.ones(len(binderGrainsDiameterList))*2))
            
                print(materialList)

                # Set minimum particle to be smallest of the two materials 
                minPar = min(grainAggMin,grainBinderMin)




        # Create empty lists if not cementStructure
        PoresDiameterList, ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0






        if multiMatToggle == "Off":

            if numCPU > 1:
            
                
                for increment in range(numIncrements-1):

                    process_pool = multiprocessing.Pool(numCPU)

                    outputMPI = process_pool.map(functools.partial(gen_particleMPI, surfaceNodes,maxParNum, minC, maxC, meshVertices, \
                        meshTets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                        maxPar,parOffset,verbose,parDiameterList,maxEdgeLength,max_dist,internalNodes), parDiameterList[particlesPlaced:particlesPlaced+math.floor(len(parDiameterList)/numIncrements)])

                    nodeMPI = np.array(outputMPI)[:,0:3]
                    diameter = np.array(outputMPI)[:,3]
                    newMaxIter = int(max(np.array(outputMPI)[:,4]))
                    maxAttempts = int(max(np.array(outputMPI)[:,5]))

                    particlesPlaced = particlesPlaced+len(np.array(outputMPI)[:,0:3])        

                    for x in range(len(nodeMPI)):

                        # Store placed particles from this increment
                        internalNodes[particlesPlaced+x,:] = nodeMPI[x,:]

                        # Obtain extents for floating bin for node to test
                        binMin = np.array(([nodeMPI[x,0]-diameter[x]/2-maxPar/2-parOffset,\
                            nodeMPI[x,1]-diameter[x]/2-maxPar/2-parOffset,nodeMPI[x,2]-\
                            diameter[x]/2-maxPar/2-parOffset]))
                        binMax = np.array(([nodeMPI[x,0]+diameter[x]/2+maxPar/2+parOffset,\
                            nodeMPI[x,1]+diameter[x]/2+maxPar/2+parOffset,nodeMPI[x,2]+\
                            diameter[x]/2+maxPar/2+parOffset]))

                        # Check if particle overlapping any just added particles (ignore first one placed)
                        if x > 0:

                            overlap = check_particleOverlapMPI(nodeMPI[x,:],diameter[x],binMin,\
                                binMax,minPar,parOffset,nodeMPI[0:x],diameter[0:x])

                            if overlap == True:

                                [newMaxIter,node,iterReq] = gen_particle(surfaceNodes,\
                                    parDiameterList[particlesPlaced+x], meshVertices, \
                                    meshTets,newMaxIter,maxIter,minPar,\
                                    maxPar,parOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes)
                                
                                internalNodes[particlesPlaced+x,:] = node[0,:]


                    self.form[5].progressBar.setValue(95*((x)/len(parDiameterList))+6) 
                    self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')


            # Generate particles for length of needed aggregate (not placed via MPI)
            for x in range(particlesPlaced,len(parDiameterList)):

                # Generate particle
                [newMaxIter,node,iterReq] = gen_particle(surfaceNodes,parDiameterList[x],meshVertices,meshTets,newMaxIter,maxIter,minPar,maxPar,\
                    parOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes)

                # Update progress bar every 1% of placement
                if x % np.rint(len(parDiameterList)/100) == 0:
                    self.form[5].progressBar.setValue(80*((x)/len(parDiameterList))+6) 

                if len(parDiameterList)<=1000:
                    # Update number particles placed every 1%
                    if x % np.rint(len(parDiameterList)/100) == 0:
                        self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')
                elif len(parDiameterList)<=10000:
                    # Update number particles placed every 0.1%
                    if x % np.rint(len(parDiameterList)/1000) == 0:
                        self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')
                else:
                    # Update number particles placed every 0.01%
                    if x % np.rint(len(parDiameterList)/10000) == 0:
                        self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')

                internalNodes[x,:] = node

            self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(len(parDiameterList)) + '/' + str(len(parDiameterList)) + ')')


            materialList = np.ones(len(parDiameterList))

            # Create empty lists if not multi-material or cementStructure
            aggGrainsDiameterList, itzGrainsDiameterList, binderGrainsDiameterList, PoresDiameterList,\
                ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0,0,0,0
            particleID = np.zeros(len(parDiameterList))




        placementTime = round(time.time() - start_time,2)   
        nParticles = len(parDiameterList) 




        ########################## End Placing Particles ##############################




    if fiberToggle in ['on','On','Y','y','Yes','yes']:
        
        fiberStartTime = time.time()

        # Use fiber data from CT scan if file exists
        if fiberFile != (0 or None or [] or ""):


            # CTScanfiber data arrangmentt
            CTScanFiberData = read_ctScan_file(fiberFile)
            CTScanFiberData = np.array(CTScanFiberData).reshape(-1,10)
            p1Fibers = CTScanFiberData[:,0:3]
            p2Fibers = CTScanFiberData[:,3:6]
            orienFibers = CTScanFiberData[:,6:9]
            fiberLengths = CTScanFiberData[:,9:10]


        # Generate fibers if no CT data
        else:


            if fiberPref<0 or fiberPref>1:

                self.form[5].statusWindow.setText('Fiber orientation strength is out of range, use 0-1')

            # Calculate number of fibers needed 
            nFiber = int(round(4*tetVolume*fiberVol/(math.pi*fiberDiameter**2*fiberLength)))

            # Initialize empty fiber nodes list outside geometry
            p1Fibers = (np.zeros((nFiber,3))+2)*maxC
            p2Fibers = (np.zeros((nFiber,3))+2)*maxC
            orienFibers = (np.zeros((nFiber,3))+2)*maxC
            fiberLengths = (np.zeros((nFiber,1)))

            # Generate fibers for number required
            for x in range(0,nFiber):
                
                if x % 100 == 0:

                    self.form[5].statusWindow.setText(str(nFiber-x) + ' Fibers Remaining')


                # Generate fiber
                [p1Fiber, p2Fiber, orienFiber, lFiber] = gen_LDPMCSL_fibers(meshVertices,meshTets,coord1,\
                    coord2,coord3,coord4,maxIter,fiberLength,maxC,maxPar,\
                    np.array([fiberOrientation1, fiberOrientation2, fiberOrientation3]),fiberPref,surfaceFaces,\
                    fiberCutting)
                p1Fibers[x,:] = p1Fiber
                p2Fibers[x,:] = p2Fiber
                orienFibers[x,:] = orienFiber
                fiberLengths[x,:] = lFiber

            fiberTime = round(time.time() - fiberStartTime,2)

            self.form[5].statusWindow.setText(str(nFiber) + ' fibers placed in ' + str(fiberTime) + ' seconds')






    ########################## Begin Tetrahedralization and Tesselation ##############################


    tetTessTimeStart = time.time()

    # Generate tetrahedralization
    self.form[5].statusWindow.setText("Status: Forming tetrahedralization.") 
    tetGen = gen_LDPMCSL_tetrahedralization(internalNodes,surfaceNodes,\
        surfaceFaces,geoName,tempPath)
    self.form[5].progressBar.setValue(89) 


    # Read in tetrahedralization
    [allNodes,allTets,allEdges] = read_LDPMCSL_tetgen(Path(tempPath + geoName \
    + '.node'),Path(tempPath + geoName + '.ele'),Path(tempPath + geoName + '.edge'))
    self.form[5].progressBar.setValue(90) 


    # Generate tesselation
    self.form[5].statusWindow.setText("Status: Forming tesselation.") 
    

    [tetFacets,facetCenters,facetAreas,facetNormals,tetn1,tetn2,tetPoints,allDiameters,facetPointData,facetCellData] = \
        gen_LDPMCSL_tesselation(allNodes,allTets,parDiameterList,minPar,geoName)    

    # If edge elements are turned on, perform edge computations
    if htcToggle in ['on','On']:
        edgeData = gen_LDPMCSL_flowEdges(htcLength,allNodes,allTets,tetPoints,maxPar,\
            meshVertices,meshTets,coord1,coord2,coord3,coord4,maxC)

    else:
        edgeData = 0




    self.form[5].progressBar.setValue(95) 
    tetTessTime = round(time.time() - tetTessTimeStart,2)   




    # Store values for unused features
    edgeMaterialList = 0
    cementStructure = 'Off'




    writeTimeStart = time.time()



    

    if multiMatToggle in ['on','On','Y','y','Yes','yes']:

        # Read in multi-material file
        [multiMatX,multiMatY,multiMatZ,multiMatRes,multiMatVoxels] = read_multiMat_file(multiMatFile)


        # Organize and store voxels of each material
        [aggVoxels,itzVoxels,binderVoxels,aggVoxelIDs] = sort_multiMat_voxels(multiMatVoxels)


    # Extend material lists for edge nodes
    particleID = np.concatenate((0*np.ones([len(allNodes)-\
        len(particleID),]),particleID))

    if multiMatToggle in ['on','On','Y','y','Yes','yes']:

        materialList = np.concatenate((2*np.ones([len(allNodes)-\
            len(materialList),]),materialList))

    elif cementStructure in ['on','On','Y','y','Yes','yes']:
        materialList = np.concatenate((edgeMaterialList,materialList))

    else:
        materialList = np.concatenate((0*np.ones([len(allNodes)-\
            len(materialList),]),materialList))    


    # For surface particles, find the nearest voxel and assign the material
    if multiMatToggle in ['on','On','Y','y','Yes','yes']:
        materialList = gen_multiMat_assign(allNodes,materialList,aggVoxels,itzVoxels,binderVoxels,internalNodes,multiMatX,multiMatY,multiMatZ,multiMatRes,minC)






    
    self.form[5].statusWindow.setText("Status: Generating facet data information.") 

    if elementType == "LDPM":
        [facetData,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial] = gen_LDPM_facetData(\
            allNodes,allTets,tetFacets,facetCenters,facetAreas,facetNormals,tetn1,\
            tetn2,materialList,multiMatRule,multiMatToggle,cementStructure,edgeMaterialList,facetCellData,particleID)
    elif elementType == "CSL":
        [facetData,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial] = gen_CSL_facetData(\
            allNodes,allEdges,allTets,tetFacets,facetCenters,facetAreas,facetNormals,tetn1,\
            tetn2,materialList,multiMatRule,multiMatToggle,cementStructure,edgeMaterialList,facetCellData,particleID)


    self.form[5].progressBar.setValue(98) 


    if ((fiberToggle in ['on','On','Y','y','Yes','yes']) and (fiberIntersections in ['on','On','Y','y','Yes','yes'])):

        self.form[5].statusWindow.setText('Determining fiber-facet intersections.')

        [FiberdataList,TotalIntersections,MaxInterPerFacet,TotalTet,TotalFiber,IntersectedFiber,projectedFacet]\
            = gen_LDPMCSL_facetfiberInt(p1Fibers,p2Fibers,fiberDiameter,fiberLengths,orienFibers,\
            geoName,allTets,allNodes,tetFacets,facetData,tetn1,tetn2,facetNormals,facetCenters)




    self.form[5].statusWindow.setText("Status: Writing external facet data file.") 
    # Create file of external triangle facets for plotting of cells
    #externalFacetsFile = externalFacetFile(facetData,meshVertices,surfaceFaces,geoName)





    # Initialize counter for number of facet materials switched
    matSwitched = 0


    # Calculate volume associated with each material (and adjust for high-order material rules)
    if multiMatToggle == "On":

        # Read in aggregate file
        try:
            [multiMatX,multiMatY,multiMatZ,multiMatRes,aggDistinctVoxels] = read_multiMat_file(aggFile)
        except:
            pass


        # Read in multi-material file
        [multiMatX,multiMatY,multiMatZ,multiMatRes,multiMatVoxels] = read_multiMat_file(multiMatFile)


        # Organize and store voxels of each material
        [aggVoxels,itzVoxels,binderVoxels,aggVoxelIDs] = sort_multiMat_voxels(multiMatVoxels)

        # Organize and store voxels of aggregate with the distinct ID file
        try:
            [aggVoxels,discard2,discard3,aggVoxelIDs] = sort_multiMat_voxels(aggDistinctVoxels)
        except:
            pass

        [itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct] = check_multiMat_matVol(subtetVol,facetMaterial,aggVoxels,itzVoxels,binderVoxels)
        
        if multiMatRule > 9:

            sortedData1 = sort_multiMat_mat(facetMaterial,facetVol1,facetVol2,particleMaterial,subtetVol)

            i = 0
            
            while (abs(itzVolFracSim-itzVolFracAct) > 0.02 or abs(itzVolFracSim-itzVolFracAct) == 0.00) and \
                abs(binderVolFracSim-binderVolFracAct) > 0.02 and \
                abs(aggVolFracSim-aggVolFracAct) > 0.02 and i < len(sortedData1):

                # Skip refinement for facets with same-material particles
                if sortedData1[i,3] != sortedData1[i,4]:
                    
                    # Refine material assignment based on volume fractions
                    sortedData = gen_multiMat_refine(sortedData1,
                        itzVolFracSim,binderVolFracSim,aggVolFracSim,\
                        itzVolFracAct,binderVolFracAct,aggVolFracAct,i)

                    # Recalculate and update volume fractions
                    [itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
                        aggVolFracAct] = check_multiMat_matVol(sortedData[:,5],sortedData[:,2],\
                        aggVoxels,itzVoxels,binderVoxels)     

                    sortedData1 = sortedData
                    matSwitched = matSwitched+1

                else:
                    pass

                i = i+1
                

            [dataList,facetMaterial] = gen_multiMat_reform(allTets,dataList,sortedData1)

        [itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
            aggVolFracAct] = check_multiMat_matVol(subtetVol,facetMaterial,\
            aggVoxels,itzVoxels,binderVoxels)

    if multiMatToggle == "Off":

        itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct,\
            PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
            PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,\
            matSwitched,multiMatRule = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0



    App.activeDocument().addObject('App::DocumentObjectGroup',dataFilesName)
    App.activeDocument().getObject(dataFilesName).Label = 'Data Files'

    App.activeDocument().addObject('App::DocumentObjectGroup',visualFilesName)
    App.activeDocument().getObject(visualFilesName).Label = 'Visualization Files'




    App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName))
    App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName))


    # If data files requested, generate them
    if dataFilesGen == True:

        self.form[5].statusWindow.setText("Status: Writing node data file.")

        mkData_nodes(geoName,tempPath,allNodes)

        self.form[5].statusWindow.setText("Status: Writing tet data file.")

        mkData_LDPMCSL_tets(geoName,tempPath,allTets)

        if elementType == "CSL":
            self.form[5].statusWindow.setText("Status: Writing edge data file.")
            mkData_LDPMCSL_edges(geoName,tempPath,allEdges)


        self.form[5].statusWindow.setText("Status: Writing facet data file.")

        # If data files requested, generate Facet File
        mkData_LDPMCSL_facets(geoName,tempPath,facetData)
        mkData_LDPMCSL_facetsVertices(geoName,tempPath,tetFacets)
        mkData_LDPMCSL_faceFacets(geoName,tempPath,surfaceNodes,surfaceFaces)

        if ((fiberToggle in ['on','On','Y','y','Yes','yes']) and (fiberIntersections in ['on','On','Y','y','Yes','yes'])):
                mkData_LDPMCSL_facetfiberInt(geoName,FiberdataList,TotalIntersections,MaxInterPerFacet,tempPath)

        

        self.form[5].statusWindow.setText("Status: Writing particle data file.")


        # Create diameters list (including zero edge particle diameters)
        allDiameters = np.concatenate((np.array([0.0,]*int(len(allNodes)-len(parDiameterList))),parDiameterList))

        # If data files requested, generate Particle Data File
        mkData_particles(allNodes,allDiameters,geoName,tempPath)


        if htcToggle in ['on','On']:

            self.form[5].statusWindow.setText("Status: Writing flow edge data file.")

            # Generate edge element data file
            mkData_LDPMCSL_flowEdges(geoName,edgeData,tempPath)

    # If visuals requested, generate them
    if visFilesGen == True:

        self.form[5].statusWindow.setText("Status: Writing visualization files.")

        # If visuals requested, generate Particle VTK File (note we only want to visualize internal nodes, hence the slicing)
        mkVtk_particles(internalNodes,parDiameterList,materialList[(len(allNodes)-len(internalNodes)):len(allNodes)],geoName,tempPath)

        # If visuals requested, generate Facet VTK File
        mkVtk_LDPMCSL_facets(geoName,tempPath,tetFacets,facetMaterial)

        # If visuals requested, generate flow edge VTK File
        if htcToggle in ['on','On']:

            mkVtk_LDPMCSL_flowEdges(geoName,edgeData,tempPath)
            mkIges_LDPMCSL_flowEdges(geoName,edgeData,tempPath)

        if fiberToggle in ['on','On','Y','y','Yes','yes']:
            mkVtk_LDPMCSL_fibers(p1Fibers,p2Fibers,fiberDiameter,fiberLengths,orienFibers,geoName,tempPath)
            mkVtk_LDPMCSL_projFacets(geoName,projectedFacet,tempPath)

        if ((fiberToggle in ['on','On','Y','y','Yes','yes']) and (fiberIntersections in ['on','On','Y','y','Yes','yes'])):
            mkVtk_LDPMCSL_nonIntFibers(p1Fibers,p2Fibers,fiberDiameter,fiberLengths,orienFibers,geoName,IntersectedFiber,tempPath)    


    # If single tet/cell visuals requested, generate them
    if singleTetGen == True:
        if elementType == "LDPM":
            mkVtk_LDPM_singleTetFacets(geoName,tempPath,tetFacets)
            mkVtk_LDPM_singleTetParticles(allNodes,allTets,allDiameters,geoName,tempPath)
            mkVtk_LDPM_singleTet(allNodes,allTets,geoName,tempPath)
            mkVtk_LDPM_singleCell(allNodes,allTets,parDiameterList,tetFacets,geoName,tempPath)
            mkPy_LDPM_singleParaview(geoName, outDir, outName, tempPath)
            mkPy_LDPM_singleParaviewLabels(geoName, tempPath)
        elif elementType == "CSL":
            pass
            mkVtk_LDPM_singleEdgeFacets(geoName,tempPath,allEdges,facetData,tetFacets)
            mkVtk_LDPM_singleEdgeParticles(allNodes,allEdges,allDiameters,geoName,tempPath)
            mkVtk_LDPM_singleEdge(allNodes,allEdges,geoName,tempPath)
        else:
            pass




    # Move files to selected output directory
    print('Moving files.')


    # List all files in temp directory
    file_names = os.listdir(tempPath)
        
    # Move all files to output directory    
    for file_name in file_names:
        shutil.move(os.path.join(tempPath, file_name), Path(outDir + outName))

    try:
        os.rename(Path(outDir + outName + '/' + geoName + '-para-mesh.vtk'),Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))
    except:
        pass
    os.remove(Path(outDir + outName + '/' + geoName + '2D.mesh'))
    os.remove(Path(outDir + outName + '/' + geoName + '.node'))
    os.remove(Path(outDir + outName + '/' + geoName + '.ele'))
    os.remove(Path(outDir + outName + '/' + geoName + '.edge'))



    print("Generated files written to: " + str(Path(outDir + outName)))




    if dataFilesGen == True:
        # Set linked object for node data file
        LDPMnodesData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMnodesData")                                     # create your object
        #LDPMnodesData.ViewObject.Proxy = IconViewProviderToFile(LDPMnodesData,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName).addObject(LDPMnodesData)
        LDPMnodesData.addProperty("App::PropertyFile",'Location','Node Data File','Location of node data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-nodes.dat'))
        
        # Set linked object for tet data file
        LDPMtetsData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMtetsData")                                     # create your object
        #LDPMtetsData.ViewObject.Proxy = IconViewProviderToFile(LDPMtetsData,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName).addObject(LDPMtetsData)
        LDPMtetsData.addProperty("App::PropertyFile",'Location','Tet Data File','Location of tets data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-tets.dat'))

        # Set linked object for facet data file
        LDPMfacetsData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMfacetsData")                                     # create your object
        #LDPMfacetsData.ViewObject.Proxy = IconViewProviderToFile(LDPMfacetsData,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName).addObject(LDPMfacetsData)
        LDPMfacetsData.addProperty("App::PropertyFile",'Location','Facet Data File','Location of facet data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-facets.dat'))

        # Set linked object for face facet data file
        LDPMfaceFacetsData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMfaceFacetsData")                                     # create your object
        App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName).addObject(LDPMfaceFacetsData)
        LDPMfaceFacetsData.addProperty("App::PropertyFile",'Location','Face Facet Data File','Location of face facet data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-faceFacets.dat'))

        # Set linked object for facet vertices data file
        LDPMfacetsVerticesData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMfacetsVerticesData")                                     # create your object
        App.getDocument(App.ActiveDocument.Name).getObject(dataFilesName).addObject(LDPMfacetsVerticesData)
        LDPMfacetsVerticesData.addProperty("App::PropertyFile",'Location','Facet Vertices Data File','Location of facet vertices data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-facetsVertices.dat'))




    if visFilesGen == True:
        # Set linked object for particle VTK file
        LDPMparticlesVTK = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMparticlesVTK")                                     # create your object
        #LDPMparticlesVTK.ViewObject.Proxy = IconViewProviderToFile(LDPMparticlesVTK,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMparticlesVTK)
        LDPMparticlesVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-particles.000.vtk'))


    if geoType in ['Cylinder', 'Cone', 'Sphere','Import CAD or Mesh']:

        if visFilesGen == True:
            # Insert mesh visualization and link mesh VTK file
            Fem.insert(str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk')),App.ActiveDocument.Name)        
            LDPMmeshVTK = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000')
            LDPMmeshVTK.Label = 'LDPMmeshVTK' 
            App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMmeshVTK)
            LDPMmeshVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))

            Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').ShapeColor = (0.80,0.80,0.80)
            Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').BackfaceCulling = False
            Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').MaxFacesShowInner = 0
            Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_mesh_000').DisplayMode = u"Faces"

        
        try:
            objName = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0].Name
            Gui.getDocument(App.ActiveDocument.Name).getObject(objName).Transparency = 50
            Gui.getDocument(App.ActiveDocument.Name).getObject(objName).ShapeColor = (0.80,0.80,0.80)
        except:
            pass

    else:
        if visFilesGen == True:
            # Set linked object for mesh VTK file
            LDPMmeshVTK = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMmeshVTK")                                     # create your object
            #LDPMmeshVTK.ViewObject.Proxy = IconViewProviderToFile(LDPMmeshVTK,os.path.join(ICONPATH,'FEMMeshICON.svg'))
            App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMmeshVTK)
            LDPMmeshVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))


        objName = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0].Name
        try:
            Gui.getDocument(App.ActiveDocument.Name).getObject(objName).Transparency = 0
            Gui.getDocument(App.ActiveDocument.Name).getObject(objName).ShapeColor = (0.80,0.80,0.80)
        except:
            pass


    if visFilesGen == True:
        # Insert facet visualization and link facet VTK file
        Fem.insert(str(Path(outDir + outName + '/' + geoName + '-para-facets.000.vtk')),App.ActiveDocument.Name)        
        LDPMfacetsVTK = App.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000')
        LDPMfacetsVTK.Label = 'LDPMfacetsVTK' 
        App.getDocument(App.ActiveDocument.Name).getObject(visualFilesName).addObject(LDPMfacetsVTK)
        LDPMfacetsVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-facets.000.vtk'))


        # Set visualization properties for facets
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').DisplayMode = u"Wireframe"
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').MaxFacesShowInner = 0
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').BackfaceCulling = False
        Gui.getDocument(App.ActiveDocument.Name).getObject(geoName + '_para_facets_000').ShapeColor = (0.36,0.36,0.36)


    # Remove initial mesh and insert final mesh visualization
    App.getDocument(App.ActiveDocument.Name).removeObject(meshName)
    meshFile = str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))
    Fem.insert(meshFile,App.ActiveDocument.Name)

    filename = os.path.basename(meshFile)
    filename, file_extension = os.path.splitext(filename)
    filename = re.sub("\.", "_", filename)
    filename = re.sub("/.", "_", filename)
    filename = re.sub("-", "_", filename)
    # If filename starts with a number, resub it with an underscore
    filename = re.sub("^\d", "_", filename)
    newMesh = App.getDocument(App.ActiveDocument.Name).getObject(filename)
    newMesh.Label = meshName




    # Set visualization properties for particle centers
    # Need to load differently for generated vs loaded meshes
    if geoType == "Import CAD or Mesh":
        filename = os.path.basename(cadFile)
        filename, file_extension = os.path.splitext(filename)
        filename = re.sub("\.", "_", filename)
        filename = re.sub("/.", "_", filename)
        filename = re.sub("-", "_", filename)
        # If filename starts with a number, resub it with an underscore
        filename = re.sub("^\d", "_", filename)
        if periodicToggle == "On":
            filename = filename[:-3] + "001"
        geoObj = App.getDocument(App.ActiveDocument.Name).getObject(filename)
        Gui.getDocument(App.ActiveDocument.Name).getObject(filename).BackfaceCulling = False
        Gui.getDocument(App.ActiveDocument.Name).getObject(filename).Transparency = 0
        Gui.getDocument(App.ActiveDocument.Name).getObject(filename).DisplayMode = u"Faces, Wireframe & Nodes"
        Gui.getDocument(App.ActiveDocument.Name).getObject(filename).ShapeColor = (0.80,0.80,0.80)
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject(filename))
        App.getDocument(App.ActiveDocument.Name).getObject(filename).Label = meshName
        GeoObject = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(filename)[0]
        GeoObject.Label = elementType + "geo" + str(0).zfill(3)
     
    else:
        Gui.getDocument(App.ActiveDocument.Name).getObject(filename).DisplayMode = u"Nodes"
        Gui.getDocument(App.ActiveDocument.Name).getObject(filename).PointSize = 3.00
        Gui.getDocument(App.ActiveDocument.Name).getObject(filename).PointColor = (0.00,0.00,0.00)









    # Set material properties, descriptions, and values based on selected constitutive equation set and material parameter set
    [materialProps, materialPropDesc, materialPropsVal] = gen_LDPMCSL_properties(constitutiveEQ, matParaSet)




    simProps = [\
        "TotalTime",\
        "TimestepSize",\
        "NumberOfThreads",\
        "NumberOutputSteps",\
        ]




    simPropDesc = [\
        "Description coming soon...",\
        "Description coming soon...",\
        "Description coming soon...",\
        "Description coming soon...",\
        ]





    # Remove unused material properties
    #App.getDocument(App.ActiveDocument.Name).getObject(materialName).removeProperty("References")

    # Add appropriate constitutive equation set
    App.getDocument(App.ActiveDocument.Name).getObject(materialName).addProperty("App::PropertyString",'ConstitutiveEquationSet','Base','Set of constitutive equations.').ConstitutiveEquationSet=constitutiveEQ



    # Add appropriate material properties
    for x in range(len(materialProps)):
        App.getDocument(App.ActiveDocument.Name).getObject(materialName).addProperty("App::PropertyString",materialProps[x],elementType+" Parameters",materialPropDesc[x])
        setattr(App.getDocument(App.ActiveDocument.Name).getObject(materialName),materialProps[x],str(materialPropsVal[x]))


    # Add appropriate simulation properties
    for x in range(len(simProps)):
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyFloat",simProps[x],"Simulation",simPropDesc[x])#.Density=0.25
    App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyEnumeration","Solver","Simulation","Solver software").Solver=['Project Chrono']
    App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyEnumeration","IntegrationScheme","Simulation","Integrator type").IntegrationScheme=['Explicit']












    self.form[5].progressBar.setValue(100) 


    if multiMatToggle == "Off":
        # Display sieve curve data
        mkDisp_sieveCurves(volFracPar, tetVolume, minPar, maxPar,fullerCoef,sieveCurveDiameter,sieveCurvePassing,parDiameterList)

    # Switch back to model window
    mw=Gui.getMainWindow()
    mdi=mw.findChild(QtGui.QMdiArea)
    mdi.activatePreviousSubWindow()

    # Switch to FEM GUI
    App.ActiveDocument.recompute()



    Gui.Control.closeDialog()
    Gui.activateWorkbench("FemWorkbench")
    FemGui.setActiveAnalysis(App.activeDocument().getObject(analysisName))
    





