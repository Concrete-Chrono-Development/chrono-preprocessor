import os
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

from PySide import QtCore, QtGui

from freecad.chronoConcrete import ICONPATH
from freecad.chronoConcrete import GUIPATH

from freecad.chronoConcrete.gui.readInputsLDPM          import readInputs
from freecad.chronoConcrete.gui.ccloadUIfile            import ccloadUIfile
from freecad.chronoConcrete.gui.ccloadUIicon            import ccloadUIicon

from freecad.chronoConcrete.generation.genAnalysis      import genAnalysis
from freecad.chronoConcrete.generation.genGeometry      import genGeometry
from freecad.chronoConcrete.generation.genSurfaceMesh   import genSurfMesh
from freecad.chronoConcrete.generation.particleVol      import particleVol
from freecad.chronoConcrete.generation.particleList     import particleList
from freecad.chronoConcrete.generation.surfMeshSize     import surfMeshSize


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



        self.form[5].progressBar.setValue(1) 


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
        print("Generating geometry")
        genGeo = genGeometry(dimensions,geoType,geoName)
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
        [vertices,edges,faces,tets] = genSurfMesh(analysisName,geoName,meshName,minPar,maxPar)
        self.form[5].progressBar.setValue(5) 






        # Calculate required volume of particles and sieve curve data
        [parVolTotal,cdf,cdf1,kappa_i] = particleVol(wcRatio,airFrac,fullerCoef,cementC,cementDensity,densityWater,\
            vertices,tets,minPar,maxPar,sieveCurveDiameter,sieveCurvePassing)


        # Calculate list of particle diameters for placement
        [maxParNum,parDiameterList] = particleList(parVolTotal,minPar,maxPar,newSieveCurveD,\
            cdf,kappa_i,NewSet,fullerCoef)

        # Calculation of surface mesh size
        maxEdgeLength = surfMeshSize(vertices,faces)






        print([maxParNum, parDiameterList])



        feminout.importVTKResults.export(ExportObjectList,FilePath)

        # Make object to store VTK files
        vtk_object = ObjectsFem.makePostVtkResult(doc,base_result,name = "VTK Files")











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