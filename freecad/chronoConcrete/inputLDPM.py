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


        self.icon = os.path.join(ICONPATH, "ldpm.svg")
        ui_file_A = os.path.join(GUIPATH, "ldpmMeshProps.ui")
        ui_file_B = os.path.join(GUIPATH, "geometry.ui")
        ui_file_C = os.path.join(GUIPATH, "aggregate.ui")   
        ui_file_D = os.path.join(GUIPATH, "mixDesign.ui")
        ui_file_E = os.path.join(GUIPATH, "additionalParameters.ui")
        ui_file_F = os.path.join(GUIPATH, "generation.ui")
        a = Gui.PySideUic.loadUi(ui_file_A)
        b = Gui.PySideUic.loadUi(ui_file_B)
        c = Gui.PySideUic.loadUi(ui_file_C)        
        d = Gui.PySideUic.loadUi(ui_file_D)          
        e = Gui.PySideUic.loadUi(ui_file_E)        
        f = Gui.PySideUic.loadUi(ui_file_F)
        self.form = [a, b, c, d, e, f]
        self.form[0].setWindowTitle("LDPM Meshing Settings")
        self.form[1].setWindowTitle("Geometry")
        self.form[2].setWindowTitle("Aggregate")        
        self.form[3].setWindowTitle("Mix Design")
        self.form[4].setWindowTitle("Additional Parameters")
        self.form[5].setWindowTitle("Model Generation") 

        self.form[0].setWindowIcon(QtGui.QIcon.fromTheme("cancel",QtGui.QIcon(":/icons/parametric/Part_Box_Parametric.svg")))
        self.form[1].setWindowIcon(QtGui.QIcon.fromTheme("cancel",QtGui.QIcon(":/icons/parametric/Part_Box_Parametric.svg")))
        self.form[2].setWindowIcon(QtGui.QIcon.fromTheme("cancel",QtGui.QIcon(":/icons/parametric/Part_Box_Parametric.svg")))
        self.form[3].setWindowIcon(QtGui.QIcon.fromTheme("cancel",QtGui.QIcon(":/icons/parametric/Part_Box_Parametric.svg")))
        self.form[4].setWindowIcon(QtGui.QIcon.fromTheme("cancel",QtGui.QIcon(":/icons/parametric/Part_Box_Parametric.svg")))
        self.form[5].setWindowIcon(QtGui.QIcon.fromTheme("cancel",QtGui.QIcon(":/icons/parametric/Part_Box_Parametric.svg")))

        QtCore.QObject.connect(self.form[0].readFileButton, QtCore.SIGNAL("clicked()"), self.openFile)
        QtCore.QObject.connect(self.form[5].pushButton, QtCore.SIGNAL("clicked()"), self.generation)





        self.minAgg = self.form[2].minAgg.value()
        self.maxAgg = self.form[2].maxAgg.value()



        self.elementType = "LDPM"
        self.constitutiveEQ = "quasiBrittle"
        self.meshName = self.elementType + "mesh"
        self.geoName = self.elementType + "geo"



    def getStandardButtons(self):
        # only show a close button
        # def accept() in no longer needed, since there is no OK button
        return int(0)



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


        # Store document
        docGui = Gui.activeDocument()




        if geoType == "Box":

            # Create a box and name it
            geo = App.ActiveDocument.addObject("Part::Box",geoName)
            geo.Label = geoName
            geo.Height = dimensions[0]
            geo.Width = dimensions[1] 
            geo.Length = dimensions[2]


        App.ActiveDocument.recompute()

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)


    def genAnalysis(self):

        # Analysis
        analysis_object = ObjectsFem.makeAnalysis(App.ActiveDocument,self.elementType)

        # Solver 
        solver_object = ObjectsFem.makeSolverCalculixCcxTools(App.ActiveDocument, "Project Chrono")
        solver_object.GeometricalNonlinearity = 'linear'
        solver_object.ThermoMechSteadyState = True
        solver_object.MatrixSolverType = 'default'
        solver_object.IterationsControlParameterTimeUse = False
        analysis_object.addObject(solver_object)

        # Store Material
        material_object = ObjectsFem.makeMaterialSolid(App.ActiveDocument, self.constitutiveEQ)
        mat = material_object.Material
        mat['Name'] = self.constitutiveEQ
        material_object.Material = mat
        analysis_object.addObject(material_object)


    def genSurfMesh(self):


        # Set up Gmsh
        femmesh_obj = ObjectsFem.makeMeshGmsh(App.ActiveDocument, self.meshName)
        App.ActiveDocument.getObject(self.meshName).CharacteristicLengthMin = self.minAgg
        App.ActiveDocument.getObject(self.meshName).CharacteristicLengthMax = self.maxAgg
        App.ActiveDocument.ActiveObject.Part = App.ActiveDocument.getObject(self.geoName)
        App.ActiveDocument.recompute()
        App.ActiveDocument.getObject(self.meshName).adjustRelativeLinks(App.ActiveDocument.getObject(self.elementType))
        App.ActiveDocument.getObject(self.elementType).addObject(App.ActiveDocument.getObject(self.meshName))

        # Run Gmsh
        gmsh_mesh = gmsh(femmesh_obj)
        error = gmsh_mesh.create_mesh()
        print(error)
        App.ActiveDocument.recompute()



    def generation(self):

        geoType = self.form[1].PrimitiveTypeCB.currentText()
        
        dimensions = []

        if geoType == "Box":
            dimensions.append(self.form[1].boxWidth.text())
            dimensions.append(self.form[1].boxLength.text())
            dimensions.append(self.form[1].boxHeight.text())



        genGeo = self.genGeometry(dimensions,geoType,self.geoName)







        genAna = genAnalysis()

        genSuf = genSurfMesh()








        self.form[5].progressBar.setValue(10) 




        gen = ParticleGen(maxAggD,minAggD,fullerCoef,wcRatio,\
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


        

    def reject(self):
        Gui.ActiveDocument.resetEdit()


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