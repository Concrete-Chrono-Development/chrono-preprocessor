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
## Primary Authors: Matthew Troemner
## ===========================================================================
##
## This file contains the function to perform the tesselation of the mesh and
## generates the facet connectivity, facet centers, facet areas, facet normals,
## and other facet data for CSL.
##
## ===========================================================================

import numpy as np

def genTesselationCSL(allNodes,allTets,parDiameter,minPar,geoName):
    
    """
    Variable List:
    --------------------------------------------------------------------------
    ### Inputs ###
    allNodes:        An array of coordinates for all nodes in the mesh.
    allTets:         An array of indices defining the tet connectivity 
    parDiameter:     An array of diameters for each particle.
    minPar:          The minimum particle diameter.
    geoName:         The name of the geometry.
    --------------------------------------------------------------------------
    ### Outputs ###
    tetFacets:       An array of indices defining the facet connectivity 
    facetCenters:    An array of coordinates for the center of each facet in the mesh.
    facetAreas:      An array of the area of each facet in the mesh.
    facetNormals:    An array of the normal vector of each facet in the mesh.
    tetNodeA:        An array of indices for the first node of each tetrahedron in the mesh.
    tetNodeB:        An array of indices for the second node of each tetrahedron in the mesh.
    tetPoints:       An array of coordinates for the center of each tetrahedron in the mesh.
    allDiameters:    An array of diameters for each node in the mesh.
    facetPointData:  Condensed array of facet nodes
    facetCellData:   Facet connectivity for condensed facet node list
    --------------------------------------------------------------------------
    """


    # Create diameters list (including fictitious edge particle diameters)
    allDiameters = np.concatenate((np.array([1.1*minPar,]*\
        int(len(allNodes)-len(parDiameter))),parDiameter))

    # Definition of Edge Points [Coordinates1,....,Coordinates6]
    edges1 = [allTets[:,0],allTets[:,1]]
    edges2 = [allTets[:,0],allTets[:,2]]
    edges3 = [allTets[:,0],allTets[:,3]]
    edges4 = [allTets[:,1],allTets[:,2]]
    edges5 = [allTets[:,1],allTets[:,3]]
    edges6 = [allTets[:,2],allTets[:,3]]

    edges = np.concatenate((edges1,edges2,edges3,edges4,edges5,edges6),\
        axis=1).T

    edges = np.sort(edges,axis=1)

    edgeNode1 = allNodes[(edges[:,0]-1).astype(int),:]
    edgeNode2 = allNodes[(edges[:,1]-1).astype(int),:]

    nodalDistance = np.linalg.norm(edgeNode1-edgeNode2, axis=1)
    edgeDistance = (nodalDistance - (allDiameters[(edges[:,0]-1).\
        astype(int)]/2) - (allDiameters[(edges[:,1]-1).astype(int)])/2)/2


    # Make unit vector from edgeNode 1 to edgeNode2, multiply vector 
    # by sum(agg1 and edgeDistance) and add to edgeNode2
    edgePoints = (edgeNode1-edgeNode2)/np.array([nodalDistance,]*3).T*\
        (np.array([allDiameters[(edges[:,1]-1).astype(int)],]*3).T/2+\
        np.array([edgeDistance,]*3).T)+edgeNode2


    # Form Edge Point List
    edgePoints = np.concatenate(np.split(edgePoints,6),axis=1)