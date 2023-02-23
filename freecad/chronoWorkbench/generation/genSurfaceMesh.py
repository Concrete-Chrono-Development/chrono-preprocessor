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
## Author: Matthew Troemner
## ===========================================================================
##
## Generate surface mesh using Gmsh and extract the vertices, edges, faces, 
## and tetrahedra information from the mesh.
##
## ===========================================================================

import FreeCAD as App
import ObjectsFem
import numpy as np
from femmesh.gmshtools import GmshTools as gmsh


def genSurfMesh(analysisName, geoName, meshName, minPar, maxPar):

    """
    Variable List:
    --------------------------------------------------------------------------
    ### Inputs ###
    analysisName: Name of the analysis object in the FreeCAD document.
    geoName:      Name of the geometry object in the FreeCAD document.
    meshName:     Name of the mesh object to be created in the document.
    minPar:       Minimum characteristic length parameter for the mesh.
    maxPar:       Maximum characteristic length parameter for the mesh.
    --------------------------------------------------------------------------
    ### Outputs ###
    vertices:     Array of vertex coordinates (shape: (num_vertices, 3))
    edges:        Array of edge node indices (shape: (num_edges, 2))
    faces:        Array of face node indices (shape: (num_faces, 3))
    tets:         Array of tetrahedron node indices (shape: (num_tets, 4))
    --------------------------------------------------------------------------
    """
    
    # Set up Gmsh
    femmesh_obj = ObjectsFem.makeMeshGmsh(App.ActiveDocument, meshName)
    # Set minimum and maximum characteristic lengths for the mesh
    App.ActiveDocument.getObject(meshName).CharacteristicLengthMin = minPar
    App.ActiveDocument.getObject(meshName).CharacteristicLengthMax = 2 * minPar
    App.ActiveDocument.getObject(meshName).MeshSizeFromCurvature = 0
    App.ActiveDocument.getObject(meshName).ElementOrder = u"1st"
    App.ActiveDocument.getObject(meshName).Algorithm2D = u"Delaunay"
    App.ActiveDocument.getObject(meshName).Algorithm3D = u"Delaunay"
    App.ActiveDocument.getObject(meshName).ElementDimension = u"3D"
    App.ActiveDocument.getObject(meshName).CoherenceMesh = True
    # Assign the geometry object to the mesh object
    App.ActiveDocument.ActiveObject.Part = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(geoName)[0]
    App.ActiveDocument.recompute()
    # Adjust relative links
    App.ActiveDocument.getObject(meshName).adjustRelativeLinks(App.ActiveDocument.getObject(analysisName))
    App.ActiveDocument.getObject(analysisName).addObject(App.ActiveDocument.getObject(meshName))

    # Run Gmsh to create the mesh
    gmsh_mesh = gmsh(femmesh_obj)
    error = gmsh_mesh.create_mesh()
    print(error)
    App.ActiveDocument.recompute()

    # Extract vertices, edges, faces, and tetrahedra information from the mesh
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

    return vertices, edges, faces, tets
