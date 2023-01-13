import FreeCAD as App
import ObjectsFem
from femmesh.gmshtools import GmshTools as gmsh



def genSurfMesh(analysisName,geoName,meshName,minPar,maxPar):


    # Set up Gmsh
    femmesh_obj = ObjectsFem.makeMeshGmsh(App.ActiveDocument, meshName)
    App.ActiveDocument.getObject(meshName).CharacteristicLengthMin = minPar
    App.ActiveDocument.getObject(meshName).CharacteristicLengthMax = maxPar
    App.ActiveDocument.getObject(meshName).ElementOrder = u"1st"
    App.ActiveDocument.ActiveObject.Part = App.ActiveDocument.getObject(geoName)
    App.ActiveDocument.recompute()
    App.ActiveDocument.getObject(meshName).adjustRelativeLinks(App.ActiveDocument.getObject(analysisName))
    App.ActiveDocument.getObject(analysisName).addObject(App.ActiveDocument.getObject(meshName))

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
