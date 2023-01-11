import os
import FreeCADGui as Gui
from freecad.chronoConcrete import ICONPATH
from freecad.chronoConcrete import GUIPATH
import FreeCAD as App
import Part
import Part,PartGui
import Mesh
import MeshPartGui, FreeCADGui
import MeshPart
import Mesh, Part, PartGui
import MaterialEditor
from femmesh.gmshtools import GmshTools as gmsh
import ObjectsFem
import FemGui
from PySide import QtCore, QtGui

from freecad.chronoConcrete.ldpmMeshing.particleGeneration import ParticleGen


class inputLDPMwindow:
    def __init__(self):

        # Load Button Icon
        #self.icon = os.path.join(ICONPATH, "ldpm.svg")

        # Load UI's for Side Panel
        ui_file_A = os.path.join(GUIPATH, "ldpmMeshProps.ui")
        ui_file_B = os.path.join(GUIPATH, "geometry.ui")
        ui_file_C = os.path.join(GUIPATH, "ldpmParticles.ui")   
        ui_file_D = os.path.join(GUIPATH, "ldpmMixDesign.ui")
        ui_file_E = os.path.join(GUIPATH, "ldpmAdditionalPara.ui")
        ui_file_F = os.path.join(GUIPATH, "generation.ui")
        a = Gui.PySideUic.loadUi(ui_file_A)
        b = Gui.PySideUic.loadUi(ui_file_B)
        c = Gui.PySideUic.loadUi(ui_file_C)        
        d = Gui.PySideUic.loadUi(ui_file_D)          
        e = Gui.PySideUic.loadUi(ui_file_E)        
        f = Gui.PySideUic.loadUi(ui_file_F)

        # Label, Load Icons, and Initialize Panels
        self.form = [a, b, c, d, e, f]
        self.form[0].setWindowTitle("LDPM Meshing Settings")
        self.form[1].setWindowTitle("Geometry")
        self.form[2].setWindowTitle("Particles")        
        self.form[3].setWindowTitle("Mix Design")
        self.form[4].setWindowTitle("Additional Parameters")
        self.form[5].setWindowTitle("Model Generation") 

        self.form[0].setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, "FEM_MaterialMechanicalNonlinear.svg"))))
        self.form[1].setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, "PartDesign_AdditiveBox.svg"))))
        self.form[2].setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, "Arch_Material_Group.svg"))))
        self.form[3].setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, "FEM_ConstraintFlowVelocity.svg"))))
        self.form[4].setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, "FEM_CreateNodesSet.svg"))))
        self.form[5].setWindowIcon(QtGui.QIcon.fromTheme("",QtGui.QIcon(os.path.join(ICONPATH, "ldpm.svg"))))



        # Connect Buttons
        QtCore.QObject.connect(self.form[0].readFileButton, QtCore.SIGNAL("clicked()"), self.openFile)


        # Run generation
        QtCore.QObject.connect(self.form[5].pushButton, QtCore.SIGNAL("clicked()"), self.generation)



    def readInputs(self):



        # Generation Type
        elementType         = "LDPM"
        meshName            = elementType + "mesh"
        geoName             = elementType + "geo"


        # Basic Settings
        if self.form[0].inelasticQuasi.isChecked():
            constitutiveEQ  = "quasiBrittle"
        else:
            constitutiveEQ  = "elastic"
        paramLocation       = self.form[0].paramLocation.text()
        numCPU              = self.form[0].numCPUbox.value()
        numIncrements       = self.form[0].numPIncBox.value()

        # Geometry Settings
        geoType             = self.form[1].geometryType.currentText()
        dimensions = []
        if geoType == "Box":
            dimensions.append(self.form[1].boxLength.text())
            dimensions.append(self.form[1].boxWidth.text())
            dimensions.append(self.form[1].boxHeight.text())
        if geoType == "Cylinder":
            dimensions.append(self.form[1].cylinderHeight.text())
            dimensions.append(self.form[1].cylinderRadius.text())
        if geoType == "Cone":
            dimensions.append(self.form[1].coneHeight.text())
            dimensions.append(self.form[1].coneRadius1.text())
            dimensions.append(self.form[1].coneRadius2.text())
        if geoType == "Sphere":
            dimensions.append(self.form[1].sphereRadius.text())
        if geoType == "Ellipsoid":
            dimensions.append(self.form[1].ellipsoidRadius1.text())
            dimensions.append(self.form[1].ellipsoidRadius2.text())
            dimensions.append(self.form[1].ellipsoidRadius3.text())
            dimensions.append(self.form[1].ellipsoidAngle1.text())
            dimensions.append(self.form[1].ellipsoidAngle2.text())
            dimensions.append(self.form[1].ellipsoidAngle3.text())
        if geoType == "Prism":
            dimensions.append(self.form[1].prismCircumradius.text())
            dimensions.append(self.form[1].prismHeight.text())
            dimensions.append(self.form[1].prismPolygon.text())

        # Particle Settings
        minPar              = self.form[2].minPar.value()
        maxPar              = self.form[2].maxPar.value()        
        fullerCoef          = self.form[2].fullerCoef.value()  
        sieveCurveDiameter  = self.form[2].sieveDiameters.text()        
        sieveCurvePassing   = self.form[2].sievePassing.text()   

        # Mix Design
        wcRatio             = self.form[3].wcRatio.value()
        densityWater        = self.form[3].waterDensity.text()
        cementC             = self.form[3].cementContent.text()
        densityCement       = self.form[3].cementDensity.text()
        airFrac1            = self.form[3].airFrac.value()
        airFrac2            = self.form[3].airFracArb.value()
 
        # Additional Parameters
        # ... Coming Soon ...


        return elementType, meshName, geoName,\
            constitutiveEQ, paramLocation, numCPU, numIncrements,\
            geoType, dimensions,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, densityCement, airFrac1, airFrac2



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




    def genGeometry(self,dimensions,geoType,geoName):







        if geoType == "Box":

            # Create a box and name it
            geo = App.ActiveDocument.addObject("Part::Box",geoName)
            geo.Label = geoName
            geo.Height = dimensions[0]
            geo.Width = dimensions[1] 
            geo.Length = dimensions[2]


        App.ActiveDocument.recompute()



    def genAnalysis(self,elementType,constitutiveEQ):

        # Analysis
        analysis_object = ObjectsFem.makeAnalysis(App.ActiveDocument,elementType)

        # Solver 
        solver_object = ObjectsFem.makeSolverCalculixCcxTools(App.ActiveDocument, "Project Chrono")
        solver_object.GeometricalNonlinearity = 'linear'
        solver_object.ThermoMechSteadyState = True
        solver_object.MatrixSolverType = 'default'
        solver_object.IterationsControlParameterTimeUse = False
        analysis_object.addObject(solver_object)

        # Store Material
        material_object = ObjectsFem.makeMaterialSolid(App.ActiveDocument, constitutiveEQ)
        mat = material_object.Material
        mat['Name'] = constitutiveEQ
        material_object.Material = mat
        analysis_object.addObject(material_object)


    def genSurfMesh(self,elementType,geoName,meshName,minPar,maxPar):

   

        # Set up Gmsh
        femmesh_obj = ObjectsFem.makeMeshGmsh(App.ActiveDocument, meshName)
        App.ActiveDocument.getObject(meshName).CharacteristicLengthMin = minPar
        App.ActiveDocument.getObject(meshName).CharacteristicLengthMax = maxPar
        App.ActiveDocument.getObject(meshName).ElementOrder = u"1st"
        App.ActiveDocument.ActiveObject.Part = App.ActiveDocument.getObject(geoName)
        App.ActiveDocument.recompute()
        App.ActiveDocument.getObject(meshName).adjustRelativeLinks(App.ActiveDocument.getObject(elementType))
        App.ActiveDocument.getObject(elementType).addObject(App.ActiveDocument.getObject(meshName))

        # Run Gmsh
        gmsh_mesh = gmsh(femmesh_obj)
        error = gmsh_mesh.create_mesh()
        print(error)
        App.ActiveDocument.recompute()

        femmesh = App.ActiveDocument.getObject(meshName).FemMesh
        #femmesh.Nodes[1]  # the first node, for all nodes ues femmesh.Nodes
        #femmesh.Volumes[0]  # the first volume, for all volumes use femmesh.Volumes
        #femmesh.getElementNodes(149) # nodes of the first volume, for all volumes use a for loop

        for v in femmesh.Edges:
            print(v) # Note that this starts after edges so number is not 1
            print(femmesh.getElementNodes(v))

        for v in femmesh.Faces:
            print(v) # Note that this starts after edges so number is not 1
            print(femmesh.getElementNodes(v))


    def generation(self):

        # Store document
        docGui = Gui.activeDocument()

        # Make new document and set view if does not exisit
        try:
            docGui.activeView().viewAxonometric()
        except:
            App.newDocument("Unnamed")
            docGui = Gui.activeDocument()
            docGui.activeView().viewAxonometric()

        
        # Read in inputs from input panel
        [elementType, meshName, geoName,\
            constitutiveEQ, paramLocation, numCPU, numIncrements,\
            geoType, dimensions,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, densityCement, airFrac1, airFrac2] = self.readInputs()
        self.form[5].progressBar.setValue(1) 

        # Generate geometry
        print("Generating geometry")
        genGeo = self.genGeometry(dimensions,geoType,geoName)
        self.form[5].progressBar.setValue(2) 

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)


        # Generate analysis objects
        genAna = self.genAnalysis(elementType,constitutiveEQ)
        self.form[5].progressBar.setValue(3) 


        # Generate surface mesh
        print("Generating surface mesh")
        genSuf = self.genSurfMesh(elementType,geoName,meshName,minPar,maxPar)
        self.form[5].progressBar.setValue(5) 








        feminout.importVTKResults.export(ExportObjectList,FilePath)

        # Make object to store VTK files
        vtk_object = ObjectsFem.makePostVtkResult(doc,base_result,name = "VTK Files")








        gen = ParticleGen(maxPar,minPar,fullerCoef,wcRatio,\
            cementC,volFracAir,q,maxIter,geoFile,aggOffsetCoeff,densityWater,\
            densityCement,dataType,output,verbose,fibers,dFiber,lFiber,vFiber,\
            fiberFile,multiMaterial,materialFile,maxGrainD,minGrainD,\
            grainFullerCoef,maxBinderD,minBinderD,binderFullerCoef,maxITZD,minITZD,\
            ITZFullerCoef,rebar,rebarFile1,dRebar1,rebarFile2,dRebar2,rebarFile3,dRebar3,edgeElements,\
            surfaceExtLength,fiberOrientation,orientationStrength,sieveCurveDiameter,sieveCurvePassing,\
            grainSieveCurveDiameter,grainSieveCurvePassing,binderSieveCurveDiameter,\
            binderSieveCurvePassing,itzSieveCurveDiameter,itzSieveCurvePassing,materialRule,\
            cementStructure,cementmaterialFile,cementStructureResolution,numberofPhases,\
            maxPoresD,minPoresD,PoresFullerCoef,PoresSieveCurveDiameter,PoresSieveCurvePassing,\
            maxClinkerD,minClinkerD,ClinkerFullerCoef,ClinkerSieveCurveDiameter,ClinkerSieveCurvePassing,\
            maxCHD,minCHD,CHFullerCoef,CHSieveCurveDiameter,CHSieveCurvePassing,\
            maxCSH_LDD,minCSH_LDD,CSH_LDFullerCoef,CSH_LDSieveCurveDiameter,CSH_LDSieveCurvePassing,\
            maxCSH_HDD,minCSH_HDD,CSH_HDFullerCoef,CSH_HDSieveCurveDiameter,CSH_HDSieveCurvePassing,\
            Mean_CSH_LD,Mean_CSH_HD,SDev_CSH_LD,SDev_CSH_HD,Alphac,SaturatedCSHDensity,\
            caeFile,periodic,prismX,prismY,prismZ,randomField,cutFiber,outputUnits,numIncrements,numCPU)











        # Switch to FEM GUI
        Gui.Control.closeDialog()
        Gui.activateWorkbench("FemWorkbench")
        Gui.Selection.addSelection(App.activeDocument().getObject(self.elementType),self.elementType)
        FemGui.setActiveAnalysis(App.activeDocument().getObject(self.elementType))
        Gui.Selection.clearSelection()
        Gui.Selection.addSelection(App.activeDocument().getObject(self.elementType),self.meshName)


        
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
                "MenuText": "LDPM Generation",
                "ToolTip" : "Generation of a standard geometry LDPM"}

    def Activated(self):

        Gui.Control.showDialog(inputLDPMwindow())

        return

    def IsActive(self):
        """Here you can define if the command must be active or not (greyed) if certain conditions
        are met or not. This function is optional."""
        return True

Gui.addCommand("inputLDPM", inputLDPM_Class())