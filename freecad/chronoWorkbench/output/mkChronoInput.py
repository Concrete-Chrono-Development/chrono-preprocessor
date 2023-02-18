## ================================================================================
## CHRONO WORKBENCH - github.com/Concrete-Chrono-Development/chrono-preprocessor
##
## Copyright (c) 2023 
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be found
## in the LICENSE file at the top level of the distribution and at
## github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
##
## ================================================================================
## Author: Matthew Troemner
## ================================================================================
##
## Description coming soon...
##
##
## ================================================================================

from pathlib import Path

def mkChronoInput(elementType, materialProps, materialPropsValues, simProps, simPropsValues, \
    nodesFilename, tetsFilename, facetsFilename, geoName, outDir, outName):

    # Generate cpp file for Project Chrono simulation
    with open(Path(outDir + outName + '/' + geoName + '-chrono.cpp'),"w") as f:                                       

        f.write("""
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
//
// ================================================================================

#include <chrono_multicore/physics/ChSystemParallel.h>
#include <chrono_multicore/physics/Ch3DOFContainer.h>
#include <chrono_multicore/solver/ChSolverMulticore.h>
#include <chrono/fea/ChElementLDPM.h>
#include <chrono_postprocess/ChStreamOutAscii.h>
#include <iostream>
#include <fstream>


using namespace chrono;


int main(int argc, char** argv) {
        """)

        # Store solver parameters
        f.write('// Solver Parameters\n')
        for x in range(len(simProps)):    
            f.write('double ' + simProps[x] + ' = ' + str(simPropsValues[x]) + ';\n')
        f.write('\n')

        # Store material parameters
        f.write('// ' + elementType + ' Material Parameters\n')
        for x in range(len(materialProps)):
            f.write('double ' + materialProps[x] + ' = ' + str(materialPropsValues[x]) + ';\n')
        f.write('\n')

        # Store files to read
        f.write('// Reference Data Files\n')
        f.write('std::ifstream nodesFilename = "' + nodesFilename + '";\n')
        f.write('std::ifstream tetsFilename = "' + tetsFilename + '";\n')
        f.write('std::ifstream facetsFilename = "' + facetsFilename + '";\n')
        f.write('\n')


        f.write("""
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



        # Read and store facets
        // Read the facets from the file facetsFilename
        std::vector<ChVector<> > facets;
        std::ifstream facetsFile(facetsFilename);
        if (!facetsFile.is_open()) {
            GetLog() << ''\"ERROR: Could not open facets file\"'';
            return;
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
    

    // Set material properties
    double rho = 1000; // density
    double E = 1e9;    // Young's modulus
    double nu = 0.3;   // Poisson's ratio
    // Additional LDPM parameters

    mesh->SetMaterialProperties(rho, E, nu);

    // Create a Chrono system
    ChSystemParallelNSC system;

    // Create a contact surface for the mesh
    auto surface = std::make_shared<ChContactSurfaceNodeCloud>();
    surface->SetNodeCloud(mesh);
    system.AddContactSurface(surface);

    // Create a Chrono solver and set solver settings
    auto solver = std::make_shared<ChSolverMulticore>();
    system.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-4);
    solver->SetOmega(0.8);
    solver->SetSharpnessAngle(45);


    // Use the parallel processing capabilities of Chrono to run the simulation on multiple cores
    int num_threads = 4;
    ChOMPUtils::SetNumThreads(num_threads);
    system.SetParallelThreadNumber(num_threads);
    system.GetSettings()->max_threads = num_threads;

    // Use the Chrono output system to write the visualization output for Paraview
    system.SetOutputMode(ChSystem::OutputMode::PARAVIEW);
    system.SetOutputFlags(ChSystem::OutputFlags::MESH);


    // Perform the simulation
    double time_step = 1e-6;
    double end_time = 1;
    int num_output_steps = 20; // output 20 times
    double output_interval = end_time / num_output_steps; // output every 1/17 of the total time
    int output_counter = 0;    // counter for output
    for (double t = 0; t < end_time; t += time_step) {
        system.DoStepDynamics(time_step);
        std::cout << ''\"Time: \"'' << t << std::endl;

        if (std::fmod(t, output_interval) < time_step) {
            // Write output files
            char filename[100];
            sprintf(filename, ''\"output_%04d.vtu\"'', static_cast<int>(t * 1000));
            ChStreamOutAscii write_vtu(filename);
            write_vtu.SetNumNodes(mesh->GetNnodes());
            write_vtu.SetNumElements(mesh->GetNelements());
            write_vtu << ''\"UnstructuredGrid\\n\"'';
            write_vtu << ''\"ASCII\\n\"'';
            write_vtu << ''\"DATASET UNSTRUCTURED_GRID\\n\"'';
            write_vtu << ''\"POINTS \"'' << mesh->GetNnodes() << ''\" float\\n\"'';
            for (int i = 0; i < mesh->GetNnodes(); ++i) {
                auto node = mesh->GetNode(i);
                write_vtu << node->GetPos().x << ''\" \"'' << node->GetPos().y << ''\" \"'' << node->GetPos().z << ''\"\\n\"'';
            }
            write_vtu << ''\"CELLS \"'' << mesh->GetNelements() << ''\" \"'' << mesh->GetNelements() * 5 << ''\"\\n\"'';
            for (int i = 0; i < mesh->GetNelements(); ++i) {
                auto element = mesh->GetElement(i);
                write_vtu << ''\"4 \"'' << element->GetNodeN(0)->GetIndex() << ''\" \"'' << element->GetNodeN(1)->GetIndex() << ''\" \"'' << element->GetNodeN(2)->GetIndex() << ''\" \"'' << element->GetNodeN(3)->GetIndex() << ''\"\\n\"'';
            }
            write_vtu << ''\"CELL_TYPES \"'' << mesh->GetNelements() << ''\"\\n\"'';
            for (int i = 0; i < mesh->GetNelements(); ++i) {
                write_vtu << ''\"10\\n\"'';
            }
            write_vtu << ''\"POINT_DATA \"'' << mesh->GetNnodes() << ''\"\\n\"'';
            write_vtu << ''\"VECTORS displacement float\\n\"'';
            for (int i = 0; i < mesh->GetNnodes(); ++i) {
                auto node = mesh->GetNode(i);
                write_vtu << node->GetD().x << ''\" \"'' << node->GetD().y << ''\" \"'' << node->GetD().z << ''\"\\n\"'';
            }

        // Increment output counter
        output_counter++;
    }
}
        """)
