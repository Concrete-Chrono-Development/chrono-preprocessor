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
##
##
## ===========================================================================

# pyright: reportMissingImports=false
from pathlib import Path
import FreeCAD as App
from femtools import membertools

def mkChronoInput(elementType, analysisName, materialProps, materialPropsDesc, materialPropsVals, simProps, simPropsValues, \
    nodesFilename, tetsFilename, facetsFilename, geoName, outDir, outName):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - elementType:      The type of element that will be used
    - analysisName:     The name of the analysis
    - materialProps:    The names of the material properties
    - materialPropsVals:The values of the material properties
    - simProps:         The names of the simulation properties
    - simPropsValues:   The values of the simulation properties
    - nodesFilename:    The name of the file containing the nodes
    - tetsFilename:     The name of the file containing the tets
    - facetsFilename:   The name of the file containing the facets
    - geoName:          The name of the geometry
    - outDir:           The output directory
    - outName:          The output name
    --------------------------------------------------------------------------
    ### Outputs ###
    - An input file for Project Chrono
    --------------------------------------------------------------------------
    """


    doc = App.ActiveDocument
    analysis = doc.getObject("LDPManalysis")

    with open(Path(outDir + outName + '/' + geoName + '-chrono.cpp'),"w") as output_file:

        output_file.write("""
// ================================================================================
// CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
//
// Copyright (c) 2023 
// All rights reserved. 
//
// Use of the code that generated this file is governed by a BSD-style license that
// can be found in the LICENSE file at the top level of the distribution and at
// github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
//
// ================================================================================
// Chrono Input File
// ================================================================================
// 
// Chrono Workbench developed by Northwestern University
//
// ================================================================================
        \n""")


        output_file.write("#include <chrono_multicore/physics/ChSystemParallel.h>\n")
        output_file.write("#include <chrono_multicore/solver/ChSolverMulticore.h>\n")
        output_file.write("#include <chrono/fea/ChMesh.h>\n")
        output_file.write("#include <chrono/fea/ChElementShellANCoutput_file.h>\n")
        output_file.write("#include <chrono/fea/ChVisualizationFEAmesh.h>\n")
        output_file.write("#include <chrono/fea/ChElementLDPM.h>\n")
        output_file.write("#include <fstream>\n")
        output_file.write("#include <iostream>\n")
        output_file.write("\n")
        output_file.write("using namespace chrono;\n")
        output_file.write("using namespace chrono::fea;\n")
        output_file.write("\n")
        output_file.write("int main(int argc, char** argv) {\n\n")


        # Store solver parameters
        output_file.write('\n    // Solver Parameters\n')
        for x in range(len(simProps)):    
            output_file.write('    double ' + simProps[x] + ' = ' + str(simPropsValues[x]) + ';\n')
        output_file.write('\n')

        # Store material parameters
        output_file.write('    // ' + elementType + ' Material Parameters\n')
        for x in range(len(materialProps)):
            output_file.write('    double ' + materialProps[x] + ' = ' + str(materialPropsVals[x]) + ';    // ' + materialPropsDesc[x] + '\n')
        output_file.write('\n')

        # Store files to read
        output_file.write('    // Reference Data Files\n')
        output_file.write('    std::ifstream nodesFilename = "' + nodesFilename + '";\n')
        output_file.write('    std::ifstream tetsFilename = "' + tetsFilename + '";\n')
        output_file.write('    std::ifstream facetsFilename = "' + facetsFilename + '";\n')
        output_file.write('\n')


        # Create the LDPM element mesh
        output_file.write("""
    // Create a tetrahedra mesh object
    std::shared_ptr<ChMesh> mesh(new ChMesh);

    // Add nodes to the mesh
    ChVector<> pos;
    while (nodesFilename >> pos.x >> pos.y >> pos.z) {
        mesh->AddNode(std::make_shared<ChNodeFEAxyz>(pos));
    }

    // Add elements to the mesh
    int n1, n2, n3, n4;
    while (tetsFilename >> n1 >> n2 >> n3 >> n4) {
        auto node1 = std::static_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(n1));
        auto node2 = std::static_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(n2));
        auto node3 = std::static_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(n3));
        auto node4 = std::static_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(n4));
        mesh->AddElement(std::make_shared<ChElementLDPM>(node1, node2, node3, node4));
    }

       
        // Read the facets from the file facetsFilename
        std::vector<ChVector<> > facets;
        std::ifstream facetsFile(facetsFilename);
        if (!facetsFile.is_open()) {
            std::cerr << \"ERROR: Could not open facets file\";
            return 0;
        }
        int Tet, IDx, IDy, IDz, mF;
        double Vol, pArea, cx, cy, cz, px, py, pz, qx, qy, qz, sx, sy, sz;
        std::string line;
        while (std::getline(facetsFile, line)) {
            if (line[0] == ''\"//\"'')
                continue;
            std::istringstream iss(line);
            if (!(iss >> Tet >> IDx >> IDy >> IDz >> Vol >> pArea >> cx >> cy >> cz >> px >> py >> pz >> qx >> qy >> qz >> sx >> sy >> sz >> mF))
                break;
            facets.push_back(ChVector<>(Tet, IDx, IDy, IDz, Vol, pArea, cx, cy, cz, px, py, pz, qx, qy, qz, sx, sy, sz, mF));
        }
        facetsFile.close();

        \n""")
        output_file.write('    mesh->SetMaterialProperties(')

        for x in range(len(materialProps)):
            output_file.write(materialProps[x] + ', ')
        output_file.write(');\n')        




        ### Area to be fixed
        #analysis.Nodes = []
        #analysis.Elements = []
        #analysis.Materials = []


        # Get the mesh to solve
        mesh = membertools.get_mesh_to_solve(App.getDocument(App.ActiveDocument.Name).getObject(analysisName))        

        Nodes = []
        Elements = []
        Materials = []


        # Create the finite element mesh
        output_file.write("    // Create the finite element mesh\n")
        output_file.write("    std::map<int, std::shared_ptr<ChNodeFEAbase>> node_map;\n")
        for node in Nodes:
            output_file.write("    auto node_{} = std::make_shared<ChNodeFEAxyz>(ChVector<>({}, {}, {}));\n".format(node.Index, node.Point[0], node.Point[1], node.Point[2]))
            output_file.write("    node_{}->SetIndex({});\n".format(node.Index, node.Index))
            output_file.write("    mesh.AddNode(node_{});\n".format(node.Index))
            output_file.write("    node_map[{}] = node_{};\n".format(node.Index, node.Index))
        for element in Elements:
            output_file.write("    {\n")
            output_file.write("        std::vector<std::shared_ptr<ChNodeFEAbase>> element_nodes;\n")
            for node in element.Nodes:
                output_file.write("        element_nodes.push_back(node_map[{}]);\n".format(node.Index))
            output_file.write("        auto element_{} = std::make_shared<ChElementShellANCF>();\n".format(element.Index))
            output_file.write("        element_{}->SetNodes(element_nodes[0], element_nodes[1], element_nodes[2], element_nodes[3]);\n".format(element.Index))
            output_file.write("        mesh.AddElement(element_{});\n".format(element.Index))
            output_file.write("    }\n")
        output_file.write("\n")

        # Set the node positions and velocities
        output_file.write("    // Set the node positions and velocities\n")
        for node in Nodes:
            output_file.write("    {\n")
            output_file.write("        auto node_ptr = std::dynamic_pointer_cast<ChNodeFEAxyz>(node_map[{}]);\n".format(node.Index))
            output_file.write("        node_ptr->SetPos({0}, {1}, {2});\n".format(node.Point[0], node.Point[1], node.Point[2]))
            output_file.write("        node_ptr->SetPos_dt({0}, {1}, {2});\n".format(node.Velocity[0], node.Velocity[1], node.Velocity[2]))
            output_file.write("    }\n")
        output_file.write("\n")

        # Set the nodal forces and boundary conditions
        output_file.write("    // Set the nodal forces and boundary conditions\n")
        for node in Nodes:
            if node.Force is not None:
                output_file.write("    {\n")
                output_file.write("        auto node_ptr = std::dynamic_pointer_cast<ChNodeFEAxyz>(node_map[{}]);\n".format(node.Index))
                force_vec = node.Force
                output_file.write("        ChVector<> force_vec({0}, {1}, {2});\n".format(force_vec[0], force_vec[1], force_vec[2]))
                output_file.write("        node_ptr->SetForce(force_vec);\n")
                output_file.write("    }\n")
        output_file.write("\n")

        # Create a material composition and add the materials to it
        output_file.write("    // Create a material composition and add the materials to it\n")
        output_file.write("    auto composition = std::make_shared<ChMaterialShellANCF>();\n")
        for material in Materials:
            if isinstance(material, Fem.ElasticMaterial): # type: ignore
                output_file.write("    composition->AddLayer({}, {}, {});\n".format(material.Thickness, material.ElasticModulus, material.PoissonRatio))
            elif isinstance(material, Fem.InelasticLinearMaterial):# type: ignore
                output_file.write("    composition->AddLayer({}, {}, {}, {});\n".format(material.Thickness, material.ElasticModulus, material.PoissonRatio, material.YieldStress))
            elif isinstance(material, Fem.InelasticTensionCompressionMaterial):# type: ignore
                output_file.write("    composition->AddLayer({}, {}, {}, {}, {});\n".format(material.Thickness, material.ElasticModulus, material.PoissonRatio, material.YieldStress, material.CompressiveYieldStress))
            elif isinstance(material, Fem.InelasticKinematicHardeningMaterial):# type: ignore
                output_file.write("    composition->AddLayer({}, {}, {}, {}, {}, {});\n".format(material.Thickness, material.ElasticModulus, material.PoissonRatio, material.YieldStress, material.KinematicHardeningModulus, material.KinematicHardeningCoeff))
            elif isinstance(material, Fem.InelasticCombinedMaterial):# type: ignore
                output_file.write("    composition->AddLayer({}, {}, {}, {}, {}, {}, {}, {});\n".format(material.Thickness, material.ElasticModulus, material.PoissonRatio, material.YieldStress, material.CompressiveYieldStress, material.KinematicHardeningModulus, material.IsotropicHardeningModulus, material.KinematicHardeningCoeff))
        output_file.write("\n")


        # Set the material composition for all elements
        output_file.write("    // Set the material composition for all elements\n")
        for element in Elements:
            output_file.write("    {\n")
            output_file.write("        auto element_shell = std::dynamic_pointer_cast<ChElementShellANCF>(mesh.GetElement({}));\n".format(element.Index))
            output_file.write("        element_shell->AddLayer(composition);\n")
            output_file.write("    }\n")
        output_file.write("\n")




        output_file.write("""\n
    // Create a Chrono system
    ChSystemParallelNSC system;

    // Create a Chrono solver and set solver settings
    auto solver = std::make_shared<ChSolverMulticore>();
    system.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-4);
    solver->SetOmega(0.8);
    solver->SetSharpnessAngle(45);

    // Use the parallel processing capabilities of Chrono to run the simulation on multiple cores\n""") 
        output_file.write("    int num_threads = " + str(int(analysis.NumberOfThreads)) + ";\n")
        output_file.write("    ChOMPUtils::SetNumThreads(num_threads);\n")
        output_file.write("    system.SetParallelThreadNumber(num_threads);\n")
        output_file.write("    system.GetSettings()->max_threads = num_threads;\n\n")


        # Set up the time integration loop and output files for Paraview at set time increments
        output_file.write("    // Set up the time integration loop and output files for Paraview at set time increments\n")
        output_file.write("    double next_output_time = 0;\n")
        output_file.write("    double output_step = 0.1;\n")
        output_file.write("    for (double t = 0; t <= TotalTime; t += TimestepSize) {\n\n")

        # Output mesh data to a VTK file at set time increments
        output_file.write("        // Output mesh data to a VTK file at set time increments\n")
        output_file.write("        if (t == 0 || t >= next_output_time) {\n")
        output_file.write("            std::string filename = \"output_\" + std::to_string(int(t*100)) + \".vtk\";\n")
        output_file.write("            std::ofstream vtk_file(filename);\n")
        output_file.write("            vtk_file << \"# vtk DataFile Version 2.0\\n\";\n")
        output_file.write("            vtk_file << \"FEM Mesh Data\\n\";\n")
        output_file.write("            vtk_file << \"ASCII\\n\";\n")
        output_file.write("            vtk_file << \"DATASET UNSTRUCTURED_GRID\\n\";\n")
        output_file.write("            vtk_file << \"POINTS \" << mesh.GetNnodes() << \" double\\n\";\n")
        output_file.write("            for (const auto& node : mesh.GetNodes()) {\n")
        output_file.write("                vtk_file << node->GetPos().x() << \" \" << node->GetPos().y() << \" \" << node->GetPos().z() << \"\\n\";\n")
        output_file.write("            }\n")
        output_file.write("            vtk_file << \"CELLS \" << mesh.GetNelements() << \" \" << (4 * mesh.GetNelements()) << \"\\n\";\n")
        output_file.write("            for (const auto& element : mesh.GetElements()) {\n")
        output_file.write("                auto element_shell = std::dynamic_pointer_cast<ChElementShellANCF>(element);\n")
        output_file.write("                if (element_shell) {\n")
        output_file.write("                    vtk_file << \"4 \" << element_shell->GetNodeN(0)->GetIndex() << \" \" << element_shell->GetNodeN(1)->GetIndex() << \" \" << element_shell->GetNodeN(2)->GetIndex() << \" \" << element_shell->GetNodeN(3)->GetIndex() << \"\\n\";\n")
        output_file.write("        }\n")
        output_file.write("    }\n\n")

        # Output the final mesh state to a VTK file
        output_file.write("    // Output the final mesh state to a VTK file\n")
        output_file.write("    std::string filename = \"output_final.vtk\";\n")
        output_file.write("    std::ofstream vtk_file(filename);\n")
        output_file.write("    vtk_file << \"# vtk DataFile Version 2.0\\n\";\n")
        output_file.write("    vtk_file << \"FEM Mesh Data\\n\";\n")
        output_file.write("    vtk_file << \"ASCII\\n\";\n")
        output_file.write("    vtk_file << \"DATASET UNSTRUCTURED_GRID\\n\";\n")
        output_file.write("    vtk_file << \"POINTS \" << mesh.GetNnodes() << \" double\\n\";\n")
        output_file.write("    for (const auto& node : mesh.GetNodes()) {\n")
        output_file.write("        vtk_file << node->GetPos().x() << \" \" << node->GetPos().y() << \" \" << node->GetPos().z() << \"\\n\";\n")
        output_file.write("    }\n")
        output_file.write("    vtk_file << \"CELLS \" << mesh.GetNelements() << \" \" << (4 * mesh.GetNelements()) << \"\\n\";\n")
        output_file.write("    for (const auto& element : mesh.GetElements()) {\n")
        output_file.write("        auto element_shell = std::dynamic_pointer_cast<ChElementShellANCF>(element);\n")
        output_file.write("        if (element_shell) {\n")
        output_file.write("            vtk_file << \"4 \" << element_shell->GetNodeN(0)->GetIndex() << \" \" << element_shell->GetNodeN(1)->GetIndex() << \" \" << element_shell->GetNodeN(2)->GetIndex() << \" \" << element_shell->GetNodeN(3)->GetIndex() << \"\\n\";\n")
        output_file.write("        }\n")
        output_file.write("    }\n")
        output_file.write("    vtk_file << \"CELL_TYPES \" << mesh.GetNelements() << \"\\n\";\n")
        output_file.write("    for (const auto& element : mesh.GetElements()) {\n")
        output_file.write("        auto element_shell = std::dynamic_pointer_cast<ChElementShellANCF>(element);\n")
        output_file.write("        if (element_shell) {\n")
        output_file.write("            vtk_file << \"9\\n\";\n")
        output_file.write("        }\n")
        output_file.write("    }\n")
        output_file.write("    vtk_file << \"POINT_DATA \" << mesh.GetNnodes() << \"\\n\";\n")
        output_file.write("    vtk_file << \"SCALARS displacement double 3\\n\";\n")
        output_file.write("    vtk_file << \"LOOKUP_TABLE default\\n\";\n")
        output_file.write("    for (const auto& node : mesh.GetNodes()) {\n")
        output_file.write("        vtk_file << node->GetPos().x() - node->GetX0().x() << \" \" << node->GetPos().y() - node->GetX0().y() << \" \" << node->GetPos().z() - node->GetX0().z() << \"\\n\";\n")
        output_file.write("    }\n")
        output_file.write("}\n")