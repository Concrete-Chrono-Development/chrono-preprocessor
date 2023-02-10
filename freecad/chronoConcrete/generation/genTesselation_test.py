import numpy as np


def genTesselation(allNodes, allTets, parDiameter, minPar, geoName):
    # Create diameters list (including fictitious edge particle diameters)
    n_fictitious_nodes = len(allNodes) - len(parDiameter)
    allDiameters = np.concatenate((np.ones(n_fictitious_nodes) * 1.1 * minPar, parDiameter))

    # Definition of Edge Points [Coordinates1,....,Coordinates6]
    edges = np.array([
        [allTets[:, 0], allTets[:, 1]],
        [allTets[:, 0], allTets[:, 2]],
        [allTets[:, 0], allTets[:, 3]],
        [allTets[:, 1], allTets[:, 2]],
        [allTets[:, 1], allTets[:, 3]],
        [allTets[:, 2], allTets[:, 3]],
    ])

    edges = np.sort(edges.reshape(-1, 2), axis=1)
    edgeNode1 = allNodes[(edges[:, 0] - 1).astype(int), :]
    edgeNode2 = allNodes[(edges[:, 1] - 1).astype(int), :]

    nodalDistance = np.linalg.norm(edgeNode1 - edgeNode2, axis=1)
    edgeDistance = (nodalDistance - (allDiameters[(edges[:, 0] - 1).astype(int)] / 2) - (
                allDiameters[(edges[:, 1] - 1).astype(int)] / 2)) / 2

    # Make unit vector from edgeNode 1 to edgeNode2, multiply vector
    # by sum(agg1 and edgeDistance) and add to edgeNode2
    edgePoints = edgeNode2 + (edgeNode1 - edgeNode2) / nodalDistance[:, np.newaxis] * (
                (allDiameters[(edges[:, 1] - 1).astype(int)][:, np.newaxis] / 2) + (
                        edgeDistance[:, np.newaxis]))

    # Form Edge Point List
    edgePoints = edgePoints.reshape(-1, 3, 6).transpose(0, 2, 1).reshape(-1, 3)

    # Definition of Face Points faceNPoint[Coordinates]

    # Face 0 (Nodes: 1,2,3) 
    nodal_distance = np.linalg.norm(allNodes[allTets[:, [1, 2, 3]] - 1] - edgePoints[:, [15, 12, 9]], axis=2)
    offset_distance = (nodal_distance - allDiameters[allTets[:, [1, 2, 3]] - 1] / 2) / 2

    face_points = (allNodes[allTets[:, [1, 2, 3]] - 1] - edgePoints[:, [15, 12, 9]]) / np.expand_dims(nodal_distance, axis=2) * np.expand_dims(offset_distance, axis=2) + edgePoints[:, [15, 12, 9]]
    face0_point = np.mean(face_points, axis=1)

    # Face 1 (Nodes: 0,2,3)
    nodal_distance = np.linalg.norm(allNodes[allTets[:, [0, 2, 3]] - 1] - edgePoints[:, [15, 6, 3]], axis=2)
    offset_distance = (nodal_distance - allDiameters[allTets[:, [0, 2, 3]] - 1] / 2) / 2

    face_points = (allNodes[allTets[:, [0, 2, 3]] - 1] - edgePoints[:, [15, 6, 3]]) / np.expand_dims(nodal_distance, axis=2) * np.expand_dims(offset_distance, axis=2) + edgePoints[:, [15, 6, 3]]
    face1_point = np.mean(face_points, axis=1)

    # Face 2 (Nodes: 0,1,3)
    nodal_distance = np.linalg.norm(allNodes[allTets[:, [0, 1, 3]] - 1] - edgePoints[:, [12, 6, 0]], axis=2)
    offset_distance = (nodal_distance - allDiameters[allTets[:, [0, 1, 3]] - 1] / 2) / 2

    face_points = (allNodes[allTets[:, [0, 1, 3]] - 1] - edgePoints[:, [12, 6, 0]]) / np.expand_dims(nodal_distance, axis=2) * np.expand_dims(offset_distance, axis=2) + edgePoints[:, [12, 6, 0]]
    face2_point = np.mean(face_points, axis=1)

    # Face 3 (Nodes: 0,1,2)
    nodal_distance = np.linalg.norm(allNodes[allTets[:, [2, 1, 0]] - 1] - edgePoints[:, [0, 3, 9]], axis=2)
    offset_distance = (nodal_distance - allDiameters[allTets[:, [2, 1, 0]] - 1] / 2) / 2

    face_points = (allNodes[allTets[:, [2, 1, 0]] - 1] - edgePoints[:, [0, 3, 9]]) / np.expand_dims(nodal_distance, axis=2) * np.expand_dims(offset_distance, axis=2) + edgePoints[:, [0, 3, 9]]
    face3_point = face_points.mean(axis=1)



    # Definition of Tet-Points [Coordinates]
    nodes = allNodes[allTets - 1]
    diameters = allDiameters[allTets - 1]
    faces = [face0Point, face1Point, face2Point, face3Point]

    tetNodalDistance = np.linalg.norm(nodes - faces, axis=2)
    tetOffsetDistance = (tetNodalDistance - diameters / 2) / 2

    tetPoints = (nodes - faces) / np.expand_dims(tetNodalDistance, axis=2) * np.expand_dims(tetOffsetDistance, axis=2) + faces
    tetPoints = tetPoints.reshape(allTets.shape[0], 12)

    tetPoint = tetPoints.mean(axis=1)


    # Facets for each vertex in General Tet 
    facet_dict = {
        1: [tetPoints, face0Point, edgePoints[:, 9:12]],
        2: [tetPoints, face0Point, edgePoints[:, 12:15]],
        3: [tetPoints, face0Point, edgePoints[:, 15:18]],
        4: [tetPoints, face1Point, edgePoints[:, 3:6]],
        5: [tetPoints, face1Point, edgePoints[:, 6:9]],
        6: [tetPoints, face1Point, edgePoints[:, 15:18]],
        7: [tetPoints, face2Point, edgePoints[:, 0:3]],
        8: [tetPoints, face2Point, edgePoints[:, 6:9]],
        9: [tetPoints, face2Point, edgePoints[:, 12:15]],
        10: [tetPoints, face3Point, edgePoints[:, 0:3]],
        11: [tetPoints, face3Point, edgePoints[:, 3:6]],
        12: [tetPoints, face3Point, edgePoints[:, 9:12]],
    }

    facet_list = [np.concatenate(facet_dict[i], axis=1) for i in range(1, 13)]

    facetsP0 = np.concatenate((facet_list[3], facet_list[4], facet_list[6], facet_list[7], facet_list[9], facet_list[10]), axis=1)
    facetsP1 = np.concatenate((facet_list[0], facet_list[1], facet_list[6], facet_list[8], facet_list[9], facet_list[11]), axis=1)
    facetsP2 = np.concatenate((facet_list[0], facet_list[2], facet_list[3], facet_list[5], facet_list[10], facet_list[11]), axis=1)
    facetsP3 = np.concatenate((facet_list[1], facet_list[2], facet_list[4], facet_list[5], facet_list[7], facet_list[8]), axis=1)

    cellTetFacets = np.concatenate((facetsP0, facetsP1, facetsP2, facetsP3), axis=1)
    tetFacets = np.concatenate(facet_list, axis=1)
    facets = tetFacets.reshape(-1, 9)

    threes = 3 * np.ones((len(facets), 1))

    facetCenters = (facets[:, 0:3] + facets[:, 3:6] + facets[:, 6:9]) / threes

    vectorAB = tetFacets.reshape(-1, 9)[:, 3:6] - tetFacets.reshape(-1, 9)[:, 0:3]
    vectorAC = tetFacets.reshape(-1, 9)[:, 6:9] - tetFacets.reshape(-1, 9)[:, 0:3]

    facetAreas = 0.5 * np.linalg.norm(np.cross(vectorAB, vectorAC), axis=1)
    facetNormals = np.cross(vectorAB, vectorAC) / np.linalg.norm(np.cross(vectorAB, vectorAC), axis=1)[:, None]

    tetn1 = [1, 1, 2, 0, 0, 2, 0, 0, 1, 0, 0, 1]
    tetn2 = [2, 3, 3, 2, 3, 3, 1, 3, 3, 1, 2, 2]

    return tetFacets, facetCenters, facetAreas, facetNormals, \
        tetn1, tetn2, tetPoints, allDiameters