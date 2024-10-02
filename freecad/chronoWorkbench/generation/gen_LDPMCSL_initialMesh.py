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
## Developed by Northwestern University
## For U.S. Army ERDC Contract No. W9132T22C0015
## Primary Authors: Matthew Troemner
## ===========================================================================
##
## Generate initial mesh using Gmsh and extract the meshVertices,  
## and tetrahedra information from the mesh.
##
## ===========================================================================

import os
import re

import FreeCAD as App #type: ignore
import ImportGui
import Fem
import ObjectsFem #type: ignore
import numpy as np
from femmesh.gmshtools import GmshTools as gmsh #type: ignore
import femmesh.femmesh2mesh as mesh2mesh #type: ignore
import Mesh #type: ignore


def gen_LDPMCSL_initialMesh(cadFile,analysisName, geoName, meshName, minPar):

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
    meshVertices:     Array of vertex coordinates (shape: (num_meshVertices, 3))
    meshTets:         Array of tetrahedron node indices (shape: (num_meshTets, 4))
    --------------------------------------------------------------------------
    """
    


    # Check if filetype is CAD (needing meshing) or mesh (already meshed)
    fileName = cadFile.split(".")
    fileExtension = fileName[-1]




    if fileExtension in ["inp", "vtk", "vtu"]:

    # If the file is a mesh file, import the mesh

        Fem.insert(cadFile,App.ActiveDocument.Name)
        filename = os.path.basename(cadFile)
        filename, file_extension = os.path.splitext(filename)
        filename = re.sub("\.", "_", filename)
        filename = re.sub("/.", "_", filename)
        filename = re.sub("-", "_", filename)
        # If filename starts with a number, resub it with an underscore
        filename = re.sub("^\d", "_", filename)
        meshObj = App.getDocument(App.ActiveDocument.Name).getObject(filename[:-3] + "001")
        meshObj.Label = meshName  


    # If the file is a CAD file, mesh the geometry
    # Or if building the mesh from scratch, create the mesh
    else:

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

    # Get mesh and initialize lists
    femmesh = App.ActiveDocument.getObjectsByLabel(meshName)[0].FemMesh
    meshVertices = []
    meshTets = []

    # Get the vertex coordinates from the mesh  
    for v in femmesh.Nodes:
        meshVertices.append(femmesh.Nodes[v])
    meshVertices = np.asarray(meshVertices)

    # Get the tetrahedra information from the mesh
    for v in femmesh.Volumes:
        meshTets.append(femmesh.getElementNodes(v))
    meshTets = np.asarray(meshTets)
    meshTets = (meshTets).astype(int)



    out_mesh = mesh2mesh.femmesh_2_mesh(App.ActiveDocument.getObjectsByLabel(meshName)[0].FemMesh)
    out_mesh = Mesh.Mesh(out_mesh)


    # Get mesh and initialize lists
    surfaceNodes = []
    surfaceFaces = []


    # Get the vertex coordinates from the mesh  
    for v in range(len(out_mesh.getPoints(0)[0])):
        surfaceNodes.append(out_mesh.getFaces(0)[0][v])
    surfaceNodes = np.asarray(surfaceNodes)

    # Get the faces information from the mesh
    for v in range(len(out_mesh.getFaces(0)[1])):
        surfaceFaces.append(out_mesh.getFaces(0)[1][v])
    surfaceFaces = np.asarray(surfaceFaces)


    return meshVertices, meshTets, surfaceNodes, surfaceFaces
