import os
import sys
import time
import tempfile
import numpy as np
from pathlib import Path
import multiprocessing
import functools
import math

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
import Spreadsheet

from PySide import QtCore, QtGui

from freecad.chronoConcrete                                     import ICONPATH
from freecad.chronoConcrete                                     import GUIPATH
from freecad.chronoConcrete                                     import TETGENPATH

from freecad.chronoConcrete.util.ccloadUIfile                   import ccloadUIfile
from freecad.chronoConcrete.util.ccloadUIicon                   import ccloadUIicon

from freecad.chronoConcrete.generation.genAnalysis              import genAnalysis
from freecad.chronoConcrete.generation.genGeometry              import genGeometry
from freecad.chronoConcrete.generation.genSurfaceMesh           import genSurfMesh
from freecad.chronoConcrete.generation.particleVol              import particleVol
from freecad.chronoConcrete.generation.particleList             import particleList
from freecad.chronoConcrete.generation.particleFaces            import particleFaces
from freecad.chronoConcrete.generation.surfMeshSize             import surfMeshSize
from freecad.chronoConcrete.generation.surfMeshExtents          import surfMeshExtents
from freecad.chronoConcrete.generation.genParticle              import generateParticle
from freecad.chronoConcrete.generation.genParticleMPI           import generateParticleMPI
from freecad.chronoConcrete.generation.particleInsideCheck      import insideCheck
from freecad.chronoConcrete.generation.particleOverlapCheck     import overlapCheck
from freecad.chronoConcrete.generation.particleOverlapCheckMPI  import overlapCheckMPI
from freecad.chronoConcrete.generation.readTetgen               import readTetgen
from freecad.chronoConcrete.generation.genTetrahedralization    import genTetrahedralization
from freecad.chronoConcrete.generation.genTesselation           import genTesselation
from freecad.chronoConcrete.generation.genFacetData             import genFacetData

from freecad.chronoConcrete.input.readInputsLDPM                import readInputs

from freecad.chronoConcrete.output.mkVtkParticles               import mkVtkParticles
from freecad.chronoConcrete.output.mkVtkFacets                  import mkVtkFacets
from freecad.chronoConcrete.output.mkDataNodes                  import mkDataNodes
from freecad.chronoConcrete.output.mkDataTets                   import mkDataTets


# Turn off error for divide by zero and invalid operations
np.seterr(divide='ignore', invalid='ignore')



#sys.executable = 'C:/Users/mtroe/anaconda3/python.exe'
#multiprocessing.set_executable('C:/Users/mtroe/anaconda3/pythonw.exe')




class inputWindow_LDPM_CSL:
    def __init__(self):

        self.form = []

        # Load UI's for Side Panel
        self.form.append(ccloadUIfile("LDPM_CSL_modelProps.ui"))
        self.form.append(ccloadUIfile("LDPM_CSL_geometry.ui"))
        self.form.append(ccloadUIfile("LDPM_CSL_particles.ui"))        
        self.form.append(ccloadUIfile("LDPM_CSL_mixDesign.ui"))          
        self.form.append(ccloadUIfile("LDPM_CSL_additionalPara.ui"))       
        self.form.append(ccloadUIfile("LDPM_CSL_generation.ui"))

        # Label, Load Icons, and Initialize Panels
        self.form[0].setWindowTitle("Model Settings")
        self.form[1].setWindowTitle("Geometry")
        self.form[2].setWindowTitle("Particles")        
        self.form[3].setWindowTitle("Mix Design")
        self.form[4].setWindowTitle("Additional Parameters")
        self.form[5].setWindowTitle("Model Generation") 

        ccloadUIicon(self.form[0],"FEM_MaterialMechanicalNonlinear.svg")
        ccloadUIicon(self.form[1],"PartDesign_AdditiveBox.svg")
        ccloadUIicon(self.form[2],"Arch_Material_Group.svg")
        ccloadUIicon(self.form[3],"FEM_ConstraintFlowVelocity.svg")
        ccloadUIicon(self.form[4],"FEM_CreateNodesSet.svg")
        ccloadUIicon(self.form[5],"ldpm.svg")

        # Set initial output directory
        self.form[5].outputDir.setText(str(Path(App.ConfigGet('UserHomePath') + '/chronoWorkbench')))

        # Connect Open File Buttons
        QtCore.QObject.connect(self.form[0].readFileButton, QtCore.SIGNAL("clicked()"), self.openFile)
        QtCore.QObject.connect(self.form[5].readDirButton, QtCore.SIGNAL("clicked()"), self.openDir)

        # Run generation for LDPM or CSL
        QtCore.QObject.connect(self.form[5].generateLDPM, QtCore.SIGNAL("clicked()"), self.generation)
        QtCore.QObject.connect(self.form[5].generateCSL, QtCore.SIGNAL("clicked()"), self.generation)
        QtCore.QObject.connect(self.form[5].writePara, QtCore.SIGNAL("clicked()"), self.generation)



    def getStandardButtons(self):

        # Only show a close button
        # def accept() in no longer needed, since there is no OK button
        return int(QtGui.QDialogButtonBox.Close)



    def openFile(self):

        path = App.ConfigGet("UserHomePath")

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getOpenFileName(None,QString.fromLocal8Bit("Read a file parameter file"),path,             "*.para") # PyQt4
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getOpenFileName(None, "Read a file parameter file", path,             "*.para") #PySide
        #                                                                     "here the text displayed on windows" "here the filter (extension)"   
        if OpenName == "":                                                            # if the name file are not selected then Abord process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            App.Console.PrintMessage("Read "+OpenName+"\n")                           # text displayed to Report view (Menu > View > Report view checked)




    def openDir(self):

        path = App.ConfigGet('UserHomePath')

        OpenName = ""
        try:
            OpenName = QtGui.QFileDialog.getExistingDirectory(None, "Open Directory",path,QtGui.QFileDialog.Option.ShowDirsOnly) 
         
        except Exception:
            OpenName, Filter = QtGui.QFileDialog.getExistingDirectory(None, "Open Directory",path,QtGui.QFileDialog.Option.ShowDirsOnly) 

        

        if OpenName == "":                                                            # if not selected then Abort process
            App.Console.PrintMessage("Process aborted"+"\n")
        else:
            self.form[5].outputDir.setText(OpenName)

        return OpenName




    def generation(self):

        # Make output directory if does not exist
        outDir =  self.form[5].outputDir.text()
        try:
            os.mkdir(outDir)
        except:
            pass

        # Initialize code start time to measure performance
        start_time = time.time()

        # Make a temporary path location
        tempPath = tempfile.gettempdir() + "/chronoConc" + str(int(np.random.uniform(1e7,1e8))) + '/'
        os.mkdir(tempPath)



        #with open(Path(tempPath + "file.txt"), "w") as f:
        #    path = "where python > " + str(Path(tempPath + "file.txt"))
        #    os.system(path)

        #with open(Path(tempPath + "file.txt"), "r") as f:
        #    output = f.readline()

        #new_string = output.replace("python", "pythonw")

        #print(new_string)





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

        # Read in inputs from input panel
        [elementType, \
            setupFile, constitutiveEQ, matParaSet, \
            numCPU, numIncrements,maxIter,placementAlg,\
            geoType, dimensions,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
            cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, airFrac2,\
            outputDir] = readInputs(self.form)


        #if fillerC > 0:
        #    airFrac = airFrac2
        #else:
        airFrac = airFrac1
        aggOffsetCoeff = 0.2                                    # Minimum distance between particles factor 
        verbose = "On"

        self.form[5].progressBar.setValue(1) 
        self.form[5].statusWindow.setText("Status: Generating objects.") 

        geoName = elementType + "geo"
        meshName = elementType + "mesh"
        analysisName = elementType + "analysis"
        materialName = elementType + "material"

        i = 0
        while App.activeDocument().getObject(geoName) != None:
            i = i+1
            geoName = elementType + "geo" + str(i)
            meshName = elementType + "mesh" + str(i)
            analysisName = elementType + "analysis" + str(i)
            materialName = elementType + "material" + str(i)



        # Generate geometry
        self.form[5].statusWindow.setText("Status: Generating geometry.") 
        genGeo = genGeometry(dimensions,geoType,geoName)
        self.form[5].progressBar.setValue(2) 

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)
        Gui.runCommand('Std_PerspectiveCamera',1)


        # Generate analysis objects
        self.form[5].statusWindow.setText("Status: Generating analysis objects.") 
        genAna = genAnalysis(analysisName,materialName)
        self.form[5].progressBar.setValue(3) 



        # Generate surface mesh
        self.form[5].statusWindow.setText("Status: Generating surface mesh.") 
        [vertices,edges,faces,tets] = genSurfMesh(analysisName,geoName,meshName,minPar,maxPar)
        self.form[5].progressBar.setValue(5) 



        # Gets extents of geometry
        [minC,maxC] = surfMeshExtents(vertices)

        self.form[5].statusWindow.setText("Status: Calculating input data.") 
        # Calculate required volume of particles and sieve curve data
        [parVolTotal,cdf,cdf1,kappa_i] = particleVol(wcRatio,airFrac,fullerCoef,cementC,cementDensity,densityWater,\
            flyashC,silicaC,scmC,flyashDensity,silicaDensity,scmDensity,\
            vertices,tets,minPar,maxPar,sieveCurveDiameter,sieveCurvePassing)

        # Temporary to skip over sieve curve option
        newSieveCurveD = 0
        NewSet = 0

        self.form[5].statusWindow.setText("Status: Calculating list of particles.") 
        # Calculate list of particle diameters for placement
        [maxParNum,parDiameterList] = particleList(parVolTotal,minPar,maxPar,newSieveCurveD,\
            cdf,kappa_i,NewSet,fullerCoef)

        # Calculation of surface mesh size
        maxEdgeLength = surfMeshSize(vertices,faces)

        # Generates points for all external triangles
        facePoints = particleFaces(vertices,faces)

        # Basic Calcs
        aggOffset = aggOffsetCoeff*minPar

        
        # Store coordinates of tets in new format
        coord1 = vertices[tets[:,0]-1]
        coord2 = vertices[tets[:,1]-1]
        coord3 = vertices[tets[:,2]-1]
        coord4 = vertices[tets[:,3]-1]









        # Initialize empty particle nodes list outside geometry
        nodes = (np.zeros((len(parDiameterList),3))+2)*maxC




        self.form[5].statusWindow.setText('Status: Placing particles into geometry. (' + str(0) + '/' + str(len(parDiameterList)) + ')') 
        # Initialize values
        newMaxIter = 2
        particlesPlaced = 0





        #numCPU = 6
        #numIncrements = 10

        if numCPU > 1:
        
            #if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
            #    print("%s Remaining." % (len(parDiameterList)))

            for increment in range(numIncrements-1):

                process_pool = multiprocessing.Pool(numCPU)

                outputMPI = process_pool.map(functools.partial(generateParticleMPI, facePoints,maxParNum, minC, maxC, vertices, \
                    tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                    maxPar,aggOffset,verbose,parDiameterList,maxEdgeLength,nodes), parDiameterList[particlesPlaced:particlesPlaced+math.floor(len(parDiameterList)/numIncrements)])

                nodeMPI = np.array(outputMPI)[:,0:3]
                diameter = np.array(outputMPI)[:,3]
                newMaxIter = int(max(np.array(outputMPI)[:,4]))
                maxAttempts = int(max(np.array(outputMPI)[:,5]))

                particlesPlaced = particlesPlaced+len(np.array(outputMPI)[:,0:3])        

                for x in range(len(nodeMPI)):

                    # Store placed particles from this increment
                    nodes[particlesPlaced+x,:] = nodeMPI[x,:]

                    # Obtain extents for floating bin for node to test
                    binMin = np.array(([nodeMPI[x,0]-diameter[x]/2-maxPar/2-aggOffset,\
                        nodeMPI[x,1]-diameter[x]/2-maxPar/2-aggOffset,nodeMPI[x,2]-\
                        diameter[x]/2-maxPar/2-aggOffset]))
                    binMax = np.array(([nodeMPI[x,0]+diameter[x]/2+maxPar/2+aggOffset,\
                        nodeMPI[x,1]+diameter[x]/2+maxPar/2+aggOffset,nodeMPI[x,2]+\
                        diameter[x]/2+maxPar/2+aggOffset]))

                    # Check if particle overlapping any just added particles (ignore first one placed)
                    if x > 0:

                        overlap = overlapCheckMPI(nodeMPI[x,:],diameter[x],binMin,\
                            binMax,minPar,aggOffset,nodeMPI[0:x],diameter[0:x])

                        if overlap == True:

                            [newMaxIter,node,iterReq] = generateParticle(particlesPlaced+x,facePoints,\
                                parDiameterList[particlesPlaced+x],maxParNum, minC, maxC, vertices, \
                                tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                                maxPar,aggOffset,'No',parDiameterList,maxEdgeLength,nodes)
                            
                            nodes[particlesPlaced+x,:] = node[0,:]


                self.form[5].progressBar.setValue(95*((x)/len(parDiameterList))+6) 
                self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')


        # Generate particles for length of needed aggregate (not placed via MPI)
        for x in range(particlesPlaced,len(parDiameterList)):

            # Generate particle
            [newMaxIter,node,iterReq] = generateParticle(x,facePoints,\
                parDiameterList[x],maxParNum, minC, maxC, vertices, \
                tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                maxPar,aggOffset,verbose,parDiameterList,maxEdgeLength,nodes)

            # Update progress bar every 1% of placement
            if x % np.rint(len(parDiameterList)/100) == 0:
                self.form[5].progressBar.setValue(80*((x)/len(parDiameterList))+6) 

            if len(parDiameterList)<=1000:
                # Update number particles placed every 1%
                if x % np.rint(len(parDiameterList)/100) == 0:
                    self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')
            else:
                # Update number particles placed every 0.1%
                if x % np.rint(len(parDiameterList)/1000) == 0:
                    self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')

            nodes[x,:] = node[0,:]

        self.form[5].statusWindow.setText("Status: Placing particles into geometry. (" + str(len(parDiameterList)) + '/' + str(len(parDiameterList)) + ')')


        materialList = np.ones(len(parDiameterList))

        placementTime = round(time.time() - start_time,2)   
        nParticles = len(parDiameterList)

        # Create empty lists if not multi-material or cementStructure
        aggGrainsDiameterList, itzDiameterList, binderDiameterList, PoresDiameterList,\
            ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0,0,0,0










        # Add points visualization and information here





        tetTessTimeStart = time.time()


        self.form[5].statusWindow.setText("Status: Forming tetrahedralization.") 
        tetGen = genTetrahedralization(nodes,vertices,\
            faces,geoName,verbose,tempPath)
        self.form[5].progressBar.setValue(89) 



        [allNodes,allTets] = readTetgen(Path(tempPath + geoName \
        + '.node'),Path(tempPath + geoName + '.ele'))
        self.form[5].progressBar.setValue(90) 



        self.form[5].statusWindow.setText("Status: Forming tesselation.") 
        [tetFacets,facetCenters,facetAreas,facetNormals,tetn1,tetn2,tetPoints,allDiameters,facetPointData,facetCellData] = \
            genTesselation(allNodes,allTets,parDiameterList,minPar,\
            geoName)
        self.form[5].progressBar.setValue(95) 
        tetTessTime = round(time.time() - tetTessTimeStart,2)   




        # Store values for unused features
        edgeMaterialList = 0
        materialRule = 0
        multiMaterial = 'Off'
        cementStructure = 'Off'




        writeTimeStart = time.time()





     
        self.form[5].statusWindow.setText("Status: Generating facet data information.") 
        [dataList,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial] = genFacetData(\
            allNodes,allTets,tetFacets,facetCenters,facetAreas,facetNormals,tetn1,\
            tetn2,materialList,materialRule,multiMaterial,cementStructure,edgeMaterialList)
        self.form[5].progressBar.setValue(98) 







        self.form[5].statusWindow.setText("Status: Writing external facet data file.") 
        # Create file of external triangle facets for plotting of cells
        #externalFacetsFile = externalFacetFile(dataList,vertices,faces,geoName)





        # Initialize counter for number of facet materials switched
        matSwitched = 0

        itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct,\
            PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
            PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,\
            matSwitched,materialRule = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0




        App.activeDocument().addObject('App::DocumentObjectGroup','dataFiles')
        App.activeDocument().dataFiles.Label = 'Data Files'

        App.activeDocument().addObject('App::DocumentObjectGroup','visualFiles')
        App.activeDocument().visualFiles.Label = 'Visualization Files'




        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject("dataFiles"))
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addObject(App.getDocument(App.ActiveDocument.Name).getObject("visualFiles"))


        print('Writing Node Data file.')

        mkDataNodes(geoName,tempPath,allNodes)

        print('Writing Tet Data file.')

        mkDataTets(geoName,tempPath,allTets)

        print('Writing Facet Data file.')

        # If data files requested, generate Facet File
        #facetFile(geoName,dataList,allTets,dataType)

        print('Writing Particle Data file.')

        # If data files requested, generate Particle Data File
        #particleData(allNodes,allTets,parDiameterList,minPar,geoName)


        print('Writing visual files.')

        # If visuals requested, generate Particle VTK File
        mkVtkParticles(nodes,parDiameterList,materialList,geoName,tempPath)

        # If visuals requested, generate Facet VTK File
        mkVtkFacets(geoName,tetFacets,dataList,facetMaterial,multiMaterial,cementStructure,tempPath,facetPointData,facetCellData)



        writeTime = round(time.time() - writeTimeStart,2)




        # Generate Log file after run
        #mkLogFile = logFile(gmshTime,nParticles,placementTime,maxPar,\
        #    minPar,fullerCoef,wcRatio,cementC,volFracAir,q,maxIter,\
        #    geoName,aggOffset,densityWater,densityCement,allTets,dataType,\
        #    tetTessTime,writeTime,geoFile,dFiber,lFiber,vFiber,fiberFile,\
        #    multiMaterial,materialFile,maxGrainD,minGrainD,grainFullerCoef,\
        #    maxBinderD,minBinderD,binderFullerCoef,maxITZD,minITZD,ITZFullerCoef,output,fibers,\
        #    itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
        #    aggVolFracAct,sieveCurveDiameter,sieveCurvePassing,matSwitched,materialRule,\
        #    cementStructure,cementmaterialFile,maxPoresD,minPoresD,PoresFullerCoef,\
        #    PoresSieveCurveDiameter,PoresSieveCurvePassing,maxClinkerD,minClinkerD,\
        #    ClinkerFullerCoef,ClinkerSieveCurveDiameter,ClinkerSieveCurvePassing,\
        #    maxCHD,minCHD,CHFullerCoef,CHSieveCurveDiameter,CHSieveCurvePassing,\
        #    maxCSH_LDD,minCSH_LDD,CSH_LDFullerCoef,CSH_LDSieveCurveDiameter,CSH_LDSieveCurvePassing,\
        #    maxCSH_HDD,minCSH_HDD,CSH_HDFullerCoef,CSH_HDSieveCurveDiameter,CSH_HDSieveCurvePassing,\
        #    PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
        #    PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,outputUnits)





        # Move files to selected output directory
        print('Moving files.')


        outName = '/' + geoName + geoType + str(i).zfill(3)
        i = 0
        while os.path.isdir(Path(outDir + outName)):
            i = i+1
            outName = '/' + geoName + geoType + str(i).zfill(3)
            
        os.rename(Path(tempPath),Path(outDir + outName))
        os.rename(Path(outDir + outName + '/' + geoName + '-para-mesh.vtk'),Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))
        os.remove(Path(outDir + outName + '/' + geoName + '2D.mesh'))
        os.remove(Path(outDir + outName + '/' + geoName + '.node'))
        os.remove(Path(outDir + outName + '/' + geoName + '.ele'))




        # Set linked object for node data file
        LDPMnodesData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMnodesData")                                     # create your object
        LDPMnodesData.ViewObject.Proxy = IconViewProviderToFile(LDPMnodesData,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject('dataFiles').addObject(LDPMnodesData)
        LDPMnodesData.addProperty("App::PropertyFile",'Location','Node Data File','Location of node data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-nodes.dat'))
        
        # Set linked object for tet data file
        LDPMtetsData = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMtetsData")                                     # create your object
        LDPMtetsData.ViewObject.Proxy = IconViewProviderToFile(LDPMtetsData,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject('dataFiles').addObject(LDPMtetsData)
        LDPMtetsData.addProperty("App::PropertyFile",'Location','Tet Data File','Location of tets data file').Location=str(Path(outDir + outName + '/' + geoName + '-data-tets.dat'))


        # Set linked object for particle VTK file
        LDPMparticlesVTK = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMparticlesVTK")                                     # create your object
        LDPMparticlesVTK.ViewObject.Proxy = IconViewProviderToFile(LDPMparticlesVTK,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject('visualFiles').addObject(LDPMparticlesVTK)
        LDPMparticlesVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-particles.000.vtk'))

        # Set linked object for mesh VTK file
        LDPMmeshVTK = App.ActiveDocument.addObject("Part::FeaturePython", "LDPMmeshVTK")                                     # create your object
        LDPMmeshVTK.ViewObject.Proxy = IconViewProviderToFile(LDPMmeshVTK,os.path.join(ICONPATH,'FEMMeshICON.svg'))
        App.getDocument(App.ActiveDocument.Name).getObject('visualFiles').addObject(LDPMmeshVTK)
        LDPMmeshVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-mesh.000.vtk'))


        # Insert facet visualization and link facet VTK file
        Fem.insert(str(Path(outDir + outName + '/' + geoName + '-para-facet.000.vtk')),App.ActiveDocument.Name)        
        LDPMfacetsVTK = App.getDocument(App.ActiveDocument.Name).getObject('LDPMgeo_para_facet_000')
        LDPMfacetsVTK.Label = 'LDPMfacetsVTK' 
        App.getDocument(App.ActiveDocument.Name).getObject('visualFiles').addObject(LDPMfacetsVTK)
        LDPMfacetsVTK.addProperty("App::PropertyFile",'Location','Paraview VTK File','Location of Paraview VTK file').Location=str(Path(outDir + outName + '/' + geoName + '-para-facet.000.vtk'))


        # Set visualization properties for facets
        Gui.getDocument(App.ActiveDocument.Name).getObject('LDPMgeo_para_facet_000').DisplayMode = u"Wireframe"
        Gui.getDocument(App.ActiveDocument.Name).getObject('LDPMgeo_para_facet_000').MaxFacesShowInner = 0
        Gui.getDocument(App.ActiveDocument.Name).getObject('LDPMgeo_para_facet_000').BackfaceCulling = False
        Gui.getDocument(App.ActiveDocument.Name).getObject('LDPMgeo_para_facet_000').ShapeColor = (0.36,0.36,0.36)


        # Set visualization properties for particle centers
        Gui.getDocument(App.ActiveDocument.Name).getObject('LDPMmesh').DisplayMode = u"Nodes"
        Gui.getDocument(App.ActiveDocument.Name).getObject('LDPMmesh').PointSize = 3.00
        Gui.getDocument(App.ActiveDocument.Name).getObject('LDPMmesh').PointColor = (0.00,0.00,0.00)






        # Store material properties

        # Property list from:
        #--- Smith, J., & Cusatis, G. (2017). Numerical analysis of projectile penetration and perforation of plain and fiber reinforced concrete slabs. 
        #--- International Journal for Numerical and Analytical Methods in Geomechanics, 41(3), 315-337.

        if elementType == 'LDPM':
            
            materialProps = [\
                "Density",\
                "CompressiveNormalModulus",\
                "PoissonsRatio",\
                "InitialHardeningModulusRatio",\
                "DensificationRatio",\
                "TransitionalStrainRatio",\
                "ShearStrengthRatio",\
                "InitialFriction",\
                "SofteningExponent",\
                "DeviatoricStrainThresholdRatio",\
                "DeviatoricDamageParameter",\
                "AsymptoticFriction",\
                "TransitionalStress",\
                "TensileStrength",\
                "TensileCharacteristicLength",\
                "ShearSofteningModulusRatio",\
                ]




            materialPropDesc = [\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                "Description coming soon...",\
                ]


        # Remove unused material properties
        App.getDocument(App.ActiveDocument.Name).getObject(materialName).removeProperty("References")

        # Add appropriate material properties
        App.getDocument(App.ActiveDocument.Name).getObject(materialName).addProperty("App::PropertyString",'ConstitutiveEquationSet','Base','Set of constitutive equations.').ConstitutiveEquationSet=constitutiveEQ




        for x in range(len(materialProps)):
            App.getDocument(App.ActiveDocument.Name).getObject(materialName).addProperty("App::PropertyFloat",materialProps[x],elementType+" Parameters",materialPropDesc[x])#.Density=0.25


       



        

        # Hide un-needed simulation properties
        

        # Add appropriate simulation properties
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyEnumeration","Solver","Simulation","Solver software").Solver=['Project Chrono']
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyEnumeration","IntegrationScheme","Simulation","Integrator type").IntegrationScheme=['Explicit','Implicit']
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyFloat","EndTime","Simulation","Simulation duration")
        App.getDocument(App.ActiveDocument.Name).getObject(analysisName).addProperty("App::PropertyString","TimestepSize","Simulation","")











        self.form[5].progressBar.setValue(100) 






        # Switch to FEM GUI
        App.ActiveDocument.recompute()



        Gui.Control.closeDialog()
        Gui.activateWorkbench("FemWorkbench")
        FemGui.setActiveAnalysis(App.activeDocument().getObject(analysisName))
        







        
    # What to do when "Close" Button Clicked
    def reject(self):
         try:
             Gui.ActiveDocument.resetEdit()
             Gui.Control.closeDialog()
         except:
             Gui.Control.closeDialog()

        
class IconViewProviderToFile:                                       # Class ViewProvider create Property view of object
    def __init__( self, obj, icon):
        self.icone = icon
        
    def getIcon(self):                                              # GetIcon
        return self.icone        


class input_LDPM_CSL_Class():
    """My new command"""

    def GetResources(self):
        return {"Pixmap"  : os.path.join(ICONPATH, "ldpm.svg"), # the name of a svg file available in the resources
                "MenuText": "LDPM/CSL Generation",
                "ToolTip" : "Generation of an LDPM or CSL geometry"}

    def Activated(self):

        Gui.Control.showDialog(inputWindow_LDPM_CSL())

        return

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True

Gui.addCommand("mod_LDPM_CSL", input_LDPM_CSL_Class())