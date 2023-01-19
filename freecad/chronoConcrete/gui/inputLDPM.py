import os
import time
import tempfile
import numpy as np
from pathlib import Path

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

from freecad.chronoConcrete                                     import ICONPATH
from freecad.chronoConcrete                                     import GUIPATH
from freecad.chronoConcrete                                     import TETGENPATH

from freecad.chronoConcrete.gui.readInputsLDPM                  import readInputs
from freecad.chronoConcrete.gui.ccloadUIfile                    import ccloadUIfile
from freecad.chronoConcrete.gui.ccloadUIicon                    import ccloadUIicon

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

from freecad.chronoConcrete.output.mkVtkParticles               import mkVtkParticles
from freecad.chronoConcrete.output.mkVtkFacets                  import mkVtkFacets

class inputLDPMwindow:
    def __init__(self):


        # Load UI's for Side Panel
        a = ccloadUIfile("ldpmMeshProps.ui")
        b = ccloadUIfile("geometry.ui")
        c = ccloadUIfile("ldpmParticles.ui")        
        d = ccloadUIfile("ldpmMixDesign.ui")          
        e = ccloadUIfile("ldpmAdditionalPara.ui")       
        f = ccloadUIfile("generation.ui")
        self.form = [a, b, c, d, e, f]

        # Label, Load Icons, and Initialize Panels
        self.form[0].setWindowTitle("LDPM Meshing Settings")
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




        # Connect Buttons
        QtCore.QObject.connect(self.form[0].readFileButton, QtCore.SIGNAL("clicked()"), self.openFile)


        # Run generation
        QtCore.QObject.connect(self.form[5].generateLDPM, QtCore.SIGNAL("clicked()"), self.generation)
        QtCore.QObject.connect(self.form[5].generateCSL, QtCore.SIGNAL("clicked()"), self.generation)




    def getStandardButtons(self):
        # only show a close button
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
            try:                                                                      # detect error to read file
                file = open(OpenName, "r")                                            # open the file selected to read (r)  # (rb is binary)
                try:                                                                  # detect error ...
                    # here your code
                    print("here your code")
                    op = OpenName.split("/")                                          # decode the path
                    op2 = op[-1].split(".")                                           # decode the file name 
                    nomF = op2[0]                                                     # the file name are isolated

                    App.Console.PrintMessage(str(nomF)+"\n")                          # the file name are displayed

                    for ligne in file:                                                # read the file
                        X  = ligne.rstrip('\n\r') #.split()                           # decode the line
                        print(X)                                                      # print the line in report view other method 
                                                                                      # (Menu > Edit > preferences... > Output window > Redirect internal Python output (and errors) to report view checked) 
                except Exception:                                                     # if error detected to read
                    App.Console.PrintError("Error read file "+"\n")                   # detect error ... display the text in red (PrintError)
                finally:                                                              # if error detected to read ... or not error the file is closed
                    file.close()                                                      # if error detected to read ... or not error the file is closed
            except Exception:                                                         # if one error detected to read file
                App.Console.PrintError("Error in Open the file "+OpenName+"\n")       # if one error detected ... display the text in red (PrintError)







    def generation(self):

        # Initialize code start time to measure performance
        start_time = time.time()

        # Make a temporary path location
        tempPath = tempfile.gettempdir() + "/chronoConc" + str(int(np.random.uniform(1e7,1e8))) + '/'
        os.mkdir(tempPath)

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
            constitutiveEQ, paramLocation, numCPU, numIncrements,placementAlg,\
            geoType, dimensions,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
            cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, airFrac2] = readInputs(self.form)
        #if fillerC > 0:
        #    airFrac = airFrac2
        #else:
        airFrac = airFrac1
        maxIter = 50000
        aggOffsetCoeff = 0.2                                    # Minimum distance between particles factor 
        verbose = "On"

        self.form[5].progressBar.setValue(1) 
        self.form[5].statusWindow.setText("Status: Generating objects.") 

        geoName = elementType + "geo"
        meshName = elementType + "mesh"
        analysisName = elementType + "analysis"

        i = 0
        while App.activeDocument().getObject(geoName) != None:
            i = i+1
            geoName = elementType + "geo" + str(i)
            meshName = elementType + "mesh" + str(i)
            analysisName = elementType + "analysis" + str(i)



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
        genAna = genAnalysis(analysisName,constitutiveEQ)
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




        self.form[5].statusWindow.setText("Status: Placing particles into geometry.") 
        # Initialize values
        newMaxIter = 2
        particlesPlaced = 0

        if numCPU > 1:
        
            if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
                print("%s Remaining." % (len(parDiameterList)))

            for increment in range(numIncrements-1):

                process_pool = multiprocessing.Pool(numCPU)

                outputMPI = process_pool.map(functools.partial(generateParticleMPI, facePoints,maxParNum, minC, maxC, vertices, \
                    tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                    maxPar,aggOffset,verbose,parDiameterList,maxEdgeLength), parDiameterList[particlesPlaced:particlesPlaced+math.floor(len(parDiameterList)/numIncrements)])

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

                            newMaxIter = generateParticle(particlesPlaced+x,facePoints,\
                                parDiameterList[particlesPlaced+x],maxParNum, minC, maxC, vertices, \
                                tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                                maxPar,aggOffset,'No',parDiameterList,maxEdgeLength,nodes)
                            
                            nodes[particlesPlaced+x,:] = node[0,:]

                # Update progress bar every 2%
                if x % np.rint(len(parDiameterList)/100) == 0:
                    self.form[5].progressBar.setValue(95*((x)/len(parDiameterList))+6) 



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

            nodes[x,:] = node[0,:]




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
        [tetFacets,facetCenters,facetAreas,facetNormals,tetn1,tetn2,tetPoints,allDiameters] = \
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










        print('Writing Mesh Data file.')

        # If data files requested, generate Mesh File
        #meshFile(geoName,allNodes,allTets)

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
        mkVtkFacets(geoName,tetFacets,dataList,facetMaterial,multiMaterial,cementStructure,tempPath)



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







        # COPY FILES TO SELECTED DIRECTORY AND CLEANUP TEMPORARY DIRECTORY







        
        Fem.insert(str(Path(tempPath + geoName + '-para-facet.000.vtk')),App.ActiveDocument.Name)


        out_mesh = femmesh.femmesh2mesh.femmesh_2_mesh(App.ActiveDocument.LDPMgeo_para_facet_000.FemMesh)

        Mesh.show(Mesh.Mesh(out_mesh))


        Gui.getDocument(App.ActiveDocument.Name).getObject('Mesh').Lighting = u"Two side"

        App.ActiveDocument.LDPMgeo_para_facet_000.ViewObject.hide()

        Gui.getDocument(App.ActiveDocument.Name).getObject('Mesh').ShapeColor = (0.36,0.36,0.36)

        Gui.getDocument(App.ActiveDocument.Name).getObject('Mesh').LineColor = (0.00,0.00,0.00)

        Gui.getDocument(App.ActiveDocument.Name).getObject('LDPMgeo').Transparency = 20

        App.ActiveDocument.removeObject('LDPMgeo_para_facet_000')


        Gui.Selection.addSelection(App.ActiveDocument.Name,'LDPMmesh')
        Gui.runCommand('Std_ToggleVisibility',0)


        App.getDocument(App.ActiveDocument.Name).getObject('Mesh').Label = "LDPMfacets"


        #feminout.importVTKResults.export(ExportObjectList,FilePath)

        # Make object to store VTK files
        #vtk_object = ObjectsFem.makePostVtkResult(doc,base_result,name = "VTK Files")









        self.form[5].progressBar.setValue(100) 


        # Switch to FEM GUI
        Gui.Control.closeDialog()
        #Gui.activateWorkbench("FemWorkbench")
        #Gui.Selection.addSelection(App.activeDocument().getObject(self.elementType),self.elementType)
        #FemGui.setActiveAnalysis(App.activeDocument().getObject(self.elementType))
        #Gui.Selection.clearSelection()
        #Gui.Selection.addSelection(App.activeDocument().getObject(self.elementType),self.meshName)









        
    # What to do when "Close" Button Clicked
    def reject(self):
         try:
             Gui.ActiveDocument.resetEdit()
             Gui.Control.closeDialog()
         except:
             Gui.Control.closeDialog()

        
        


class inputLDPM_Class():
    """My new command"""

    def GetResources(self):
        return {"Pixmap"  : os.path.join(ICONPATH, "ldpm.svg"), # the name of a svg file available in the resources
                "MenuText": "LDPM/CSL Generation",
                "ToolTip" : "Generation of an LDPM or CSL geometry"}

    def Activated(self):

        Gui.Control.showDialog(inputLDPMwindow())

        return

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True

Gui.addCommand("inputLDPM", inputLDPM_Class())