import FreeCAD as App
import ObjectsFem



def genAnalysis(analysisName,constitutiveEQ):

    # Analysis
    analysis_object = ObjectsFem.makeAnalysis(App.ActiveDocument,analysisName)

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