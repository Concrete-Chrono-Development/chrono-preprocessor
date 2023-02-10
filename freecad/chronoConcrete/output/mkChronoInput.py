from pathlib import Path





def mkChronoInput(elementType, materialProps, materialPropsValues, simProps, simPropsValues, \
    nodesFilename, tetsFilename, geoName, outDir, outName):

    # Generate cpp file for Project Chrono simulation
    with open(Path(outDir + outName + '/' + geoName + '-Chrono.cpp'),"w") as f:                                       

        f.write('#include "chrono/physics/ChSystemSMC.h"\n')
        f.write('#include "chrono/solver/ChIterativeSolverLS.h"\n')
        f.write('#include "chrono/fea/ChElementTetraCorot_4.h"\n')
        f.write('#include "chrono/fea/ChMesh.h"\n')
        f.write('\n')
        f.write('using namespace chrono;\n')
        f.write('using namespace fea;\n')
        f.write('\n')

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
        f.write('std::string nodesFilename = "' + nodesFilename + '";\n')
        f.write('std::string tetsFilename = "' + tetsFilename + '";\n')
        f.write('\n')

        # Start main simulation section
        f.write('void mainSim() {\n')
        f.write('    GetLog() << "\\n-------------------------------------------------\\n";\n')
        f.write('    GetLog() << "Begin FEM dynamics simulation with implicit integration \\n\\n";\n')
        f.write('\n')
        f.write('    // The physical system: it contains all physical objects.\n')
        f.write('    ChSystemSMC sys;\n')
        f.write('\n')
        f.write('    // Create a mesh, that is a container for groups\n')
        f.write('    // of elements and their referenced nodes.\n')
        f.write('    auto tetMesh = chrono_types::make_shared<ChMesh>();\n')
        f.write('\n')

        # Create the material
        f.write('    // Create a material, that must be assigned to each element,\n')
        f.write('    // and set its parameters\n')
        f.write('    auto ' + elementType + 'Mat = chrono_types::make_shared<Ch' + elementType + '>();\n')
        for x in range(len(materialProps)):
            f.write('    ' + elementType + 'Mat->Set_' + materialProps[x] + '(' + materialProps[x] + ');\n')
        f.write('\n')

        # Read and store nodes
        f.write('    // Read the nodes from the file nodesFilename\n')
        f.write('    std::vector<ChVector<> > nodes;\n')
        f.write('    std::ifstream nodesFile(nodesFilename);\n')
        f.write('    if (!nodesFile.is_open()) {\n')
        f.write('        GetLog() << "ERROR: Could not open nodes file\\n";\n')
        f.write('        return;\n')
        f.write('    }\n')
        f.write('    double x, y, z;\n')
        f.write('    std::string line;\n')
        f.write('    while (std::getline(nodesFile, line)) {\n')
        f.write('        if (line[0] == ''\"//\"'')\n')
        f.write('            continue;\n')
        f.write('        std::istringstream iss(line);\n')
        f.write('        if (!(iss >> x >> y >> z))\n')
        f.write('            break;\n')
        f.write('        nodes.push_back(ChVector<>(x, y, z));\n')
        f.write('    }\n')
        f.write('    nodesFile.close()\n')
        f.write('\n')




        # Read and store tets
        f.write('    // Read the tetrahedra from the file tetsFilename\n')
        f.write('    std::vector<std::array<int, 4> > tetrahedra;\n')
        f.write('    std::ifstream tetrahedraFile(tetsFilename);\n')
        f.write('    if (!tetrahedraFile.is_open()) {\n')
        f.write('        GetLog() << "ERROR: Could not open tetrahedra file\\n";\n')
        f.write('        return;\n')
        f.write('    }\n')
        f.write('    int i, j, k, l;\n')
        f.write('    while (tetrahedraFile >> i >> j >> k >> l) {\n')
        f.write('        tetrahedra.push_back({i, j, k, l});\n')
        f.write('    }\n')
        f.write('    tetrahedraFile.close();\n')
        f.write('\n')

        # Create node objects
        f.write('    // Create ChNodeFEAxyz objects for each node and add to the mesh\n')
        f.write('    std::vector<std::shared_ptr<ChNodeFEAxyz> > mnodes;\n')
        f.write('    for (int i = 0; i < nodes.size(); i++) {\n')
        f.write('        auto mnode = chrono_types::make_shared<ChNodeFEAxyz>(nodes[i]);\n')
        f.write('        tetMesh->AddNode(mnode);\n')
        f.write('        mnodes.push_back(mnode);\n')
        f.write('    }\n')
        f.write('\n')

        # Create tet objects
        f.write('    // Create ChElementTetraCorot_4 elements for each tetrahedron\n')
        f.write('    for (int i = 0; i < tetrahedra.size(); i++) {\n')
        f.write('        auto melement = chrono_types::make_shared<ChElementTetraCorot_4>();\n')
        f.write('        melement->SetNodes(mnodes[tetrahedra[i][0]], mnodes[tetrahedra[i][1]], mnodes[tetrahedra[i][2]], mnodes[tetrahedra[i][3]]);\n')
        f.write('        melement->SetMaterial(ldpmMat);\n')
        f.write('        tetMesh->AddElement(melement);\n')
        f.write('    }\n')
        f.write('\n')

        # Add the mesh to the system
        f.write('    // Add the mesh to the system\n')
        f.write('    sys.Add(tetMesh);\n')
        f.write('\n')

        # Set up the solver
        f.write('    // Setup for a dynamic time integration:\n')
        f.write('    auto solver = chrono_types::make_shared<ChSolverMINRES>();\n')
        f.write('    sys.SetSolver(solver);\n')
        f.write('    solver->SetMaxIterations(maxIter);\n')
        f.write('    solver->SetTolerance(tolerance);\n')
        f.write('    sys.SetSolverForceTolerance(solverForceTolerance);\n')
        f.write('    sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);\n')
        f.write('\n')

        # Run the analysis
        f.write('    while (sys.GetChTime() < endTime) {\n')
        f.write('        sys.DoStepDynamics(timestep);\n')
        f.write('        GetLog() << " t =" << sys.GetChTime() << "  \\n";\n')
        f.write('    }\n')
        f.write('}\n')
        f.write('\n')

        # Form the main function
        f.write('int main(int argc, char* argv[]) {\n')
        f.write('    GetLog() << "Using Chrono version: " << CHRONO_VERSION << "\\n\\n";\n')
        f.write('\n')
        f.write('    mainSim();\n')
        f.write('    return 0;\n')
        f.write('}\n')