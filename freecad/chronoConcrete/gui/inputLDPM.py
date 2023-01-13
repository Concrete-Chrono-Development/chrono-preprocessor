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

import ObjectsFem
import FemGui
from PySide import QtCore, QtGui

from freecad.chronoConcrete.ldpmMeshing.particleGeneration import ParticleGen
from freecad.chronoConcrete.gui.readInputsLDPM import readInputs

from freecad.chronoConcrete.generation.genAnalysis import genAnalysis
from freecad.chronoConcrete.generation.genSurfaceMesh import genSurfMesh

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




    def genGeometry(self,dimensions,geoType,geoName):


        if all(float(i.strip(" mm")) > 0 for i in dimensions):
            pass
        else:
            raise Exception("One or more geometry dimensions are less than or equal to zero. Please revise.")




        if geoType == "Box":

            # Create a box and name it
            geo = App.ActiveDocument.addObject("Part::Box",geoName)
            geo.Label = geoName
            geo.Height = dimensions[0]
            geo.Width = dimensions[1]
            geo.Length = dimensions[2]


        App.ActiveDocument.recompute()









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
        Gui.runCommand('Std_PerspectiveCamera',1)

        # Read in inputs from input panel
        [elementType, \
            constitutiveEQ, paramLocation, numCPU, numIncrements,placementAlg,\
            geoType, dimensions,\
            minPar, maxPar, fullerCoef, sieveCurveDiameter, sieveCurvePassing,\
            wcRatio, densityWater, cementC, flyashC, silicaC, scmC,\
            cementDensity, flyashDensity, silicaDensity, scmDensity, airFrac1, airFrac2] = readInputs(self.form)
        self.form[5].progressBar.setValue(1) 


        geoName = elementType + "geo"
        meshName = elementType + "mesh"
        analysisName = elementType + "ana"

        i = 0
        while App.activeDocument().getObject(geoName) != None:
            i = i+1
            geoName = elementType + "geo" + str(i)
            meshName = elementType + "mesh" + str(i)
            analysisName = elementType + "ana" + str(i)



        # Generate geometry
        print("Generating geometry")
        genGeo = self.genGeometry(dimensions,geoType,geoName)
        self.form[5].progressBar.setValue(2) 

        # Set view
        docGui.activeView().viewAxonometric()
        Gui.SendMsgToActiveView("ViewFit")
        Gui.runCommand('Std_DrawStyle',6)
        Gui.runCommand('Std_PerspectiveCamera',1)

        # Generate analysis objects
        genAna = genAnalysis(analysisName,constitutiveEQ)
        self.form[5].progressBar.setValue(3) 


        # Generate surface mesh
        print("Generating surface mesh")
        genSuf = genSurfMesh(analysisName,geoName,meshName,minPar,maxPar)
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