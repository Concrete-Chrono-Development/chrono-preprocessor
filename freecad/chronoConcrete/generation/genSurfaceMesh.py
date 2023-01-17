import FreeCAD as App
import ObjectsFem
import numpy as np
from femmesh.gmshtools import GmshTools as gmsh



def genSurfMesh(analysisName,geoName,meshName,minPar,maxPar):


    # Set up Gmsh
    femmesh_obj = ObjectsFem.makeMeshGmsh(App.ActiveDocument, meshName)
    App.ActiveDocument.getObject(meshName).CharacteristicLengthMin = minPar
    App.ActiveDocument.getObject(meshName).CharacteristicLengthMax = 2*minPar
    App.ActiveDocument.getObject(meshName).MeshSizeFromCurvature = 0
    App.ActiveDocument.getObject(meshName).ElementOrder = u"1st"
    App.ActiveDocument.getObject(meshName).Algorithm2D = u"Delaunay"
    App.ActiveDocument.getObject(meshName).Algorithm3D = u"Delaunay"
    App.ActiveDocument.getObject(meshName).ElementDimension = u"3D"
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


    vertices = []
    edges = []
    faces = []
    tets = []


    for v in femmesh.Nodes:
        vertices.append(femmesh.Nodes[v])
    vertices = np.asarray(vertices)

    for v in femmesh.Edges:
        edges.append(femmesh.getElementNodes(v))
    edges = np.asarray(edges)

    for v in femmesh.Faces:
        faces.append(femmesh.getElementNodes(v))
    faces = np.asarray(faces)

    for v in femmesh.Volumes:
        tets.append(femmesh.getElementNodes(v))
    tets = np.asarray(tets)
    tets = (tets).astype(int)

    return vertices,edges,faces,tets
