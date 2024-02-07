
import numpy as np
import math
import time
import multiprocessing
import functools
from pathlib import Path


import sys
import FreeCAD as App

from freecad.chronoWorkbench.input.read_multiMat_file                     import read_multiMat_file
from check_multiMat_size                                    import check_multiMat_size
from sort_multiMat_voxels                                   import sort_multiMat_voxels
from calc_LDPMCSL_sieveCurve                                import calc_LDPMCSL_sieveCurve
from calc_LDPMCSL_parVolume                                 import calc_LDPMCSL_parVolume
from gen_LDPMCSL_particleList                               import gen_LDPMCSL_particleList
from gen_LDPMCSL_subParticle                                import gen_LDPMCSL_subParticle
from gen_LDPMCSL_particleMPI                                import gen_LDPMCSL_particleMPI
from check_LDPMCSL_particleOverlapMPI                       import check_LDPMCSL_particleOverlapMPI
from gen_LDPMCSL_particle                                   import gen_LDPMCSL_particle



def gen_multiStep(tempPath, numCPU, numIncrements, maxIter, aggOffset, maxEdgeLength, max_dist, minPar, maxPar, sieveCurveDiameter, sieveCurvePassing, wcRatio, cementC, airFrac, fullerCoef, flyashC, silicaC, scmC, fillerC, flyashDensity, silicaDensity, scmDensity, fillerDensity, cementDensity, densityWater, multiMatToggle, multiMatFile, grainAggMin, grainAggMax, grainAggFuller, grainAggSieveD, grainAggSieveP, grainBinderMin, grainBinderMax, grainBinderFuller, grainBinderSieveD, grainBinderSieveP, grainITZMin, grainITZMax, grainITZFuller, grainITZSieveD, grainITZSieveP, tetVolume, minC, maxC, verbose):

    # Load back in these seven matrices from their temporary files:
    # coord1, coord2, coord3, coord4, meshVertices, meshTets, surfaceNodes



    maxC = np.array(maxC)
    minC = np.array(minC)

    # Read in the seven matrices from their temporary files
    coord1 = np.load(tempPath + 'coord1.npy')
    coord2 = np.load(tempPath + 'coord2.npy')
    coord3 = np.load(tempPath + 'coord3.npy')
    coord4 = np.load(tempPath + 'coord4.npy')
    meshVertices = np.load(tempPath + 'meshVertices.npy')
    meshTets = np.load(tempPath + 'meshTets.npy')
    surfaceNodes = np.load(tempPath + 'surfaceNodes.npy')


    if multiMatToggle == "On":

        # Read in multi-material file
        [multiMatX,multiMatY,multiMatZ,multiMatRes,multiMatVoxels] = read_multiMat_file(multiMatFile)

        # Confirm if the voxelated multi-material file is larger than the provided geometry
        topoCheck = check_multiMat_size(multiMatX,multiMatY,multiMatZ,multiMatRes,minC,maxC)

        # Organize and store voxels of each material
        [aggVoxels,itzVoxels,binderVoxels] = sort_multiMat_voxels(multiMatVoxels)




        # Do calculations for aggregate, binder, and ITZ
        for i in range(3):
            
            if i == 0:
                [grainMin,grainMax,grainFuller,grainSieveD,grainSieveP] = [grainAggMin,grainAggMax,grainAggFuller,grainAggSieveD,grainAggSieveP]
            elif i == 1:
                [grainMin,grainMax,grainFuller,grainSieveD,grainSieveP] = [grainBinderMin,grainBinderMax,grainBinderFuller,grainBinderSieveD,grainBinderSieveP]
            elif i == 2:
                [grainMin,grainMax,grainFuller,grainSieveD,grainSieveP] = [grainITZMin,grainITZMax,grainITZFuller,grainITZSieveD,grainITZSieveP]


            # Shift sieve curve if needed
            if grainSieveD != (0 or None or [] or ""):
                [newGrainSieveCurveD,newGrainSieveCurveP,grainNewSet,grainW_min,grainW_max] = calc_LDPMCSL_sieveCurve(grainMin,grainMax,grainSieveD,grainSieveP)
            else:
                newGrainSieveCurveD,newGrainSieveCurveP,grainNewSet,grainW_min,grainW_max = 0, 0, 0, 0, 0

            # Calculates volume of each set of grains
            [volGrainFracPar,volGrains,cdf,cdf1,kappa_i] = calc_LDPMCSL_parVolume(tetVolume*len(aggVoxels)/(len(aggVoxels)+len(itzVoxels)+len(binderVoxels)), wcRatio, cementC,
                                                        airFrac, grainFuller, 
                                                        flyashC, silicaC, scmC, fillerC,
                                                        flyashDensity, silicaDensity, 
                                                        scmDensity, fillerDensity, cementDensity,
                                                        densityWater, grainMin, grainMax,
                                                        newGrainSieveCurveD, newGrainSieveCurveP, 
                                                        grainNewSet, grainW_min, grainW_max)

            # Generates list of needed grains
            [maxGrainsNum,grainsDiameterList] = gen_LDPMCSL_particleList(volGrains,grainMin,grainMax,newGrainSieveCurveD,cdf,kappa_i,grainNewSet,grainFuller)

            if i == 0:
                aggGrainsDiameterList = grainsDiameterList
            elif i == 1:
                binderGrainsDiameterList = grainsDiameterList
            elif i == 2:
                itzGrainsDiameterList = grainsDiameterList

        # Combine all grain lists (in order of aggregate > ITZ > binder -- the order they will be placed in the geometry)
        parDiameterList = np.concatenate((aggGrainsDiameterList,itzGrainsDiameterList,binderGrainsDiameterList))

        # Initialize empty list of all nodes outside geometry
        internalNodes = (np.zeros((len(aggGrainsDiameterList)+\
            len(binderGrainsDiameterList)+len(itzGrainsDiameterList),3))+2)*maxC




    if multiMatToggle == "Off":


        # Shift sieve curve if needed
        if sieveCurveDiameter != (0 or None or [] or ""):
            # Shifts sieve curve to appropriate range
            [newSieveCurveD, newSieveCurveP, NewSet, w_min, w_max] = calc_LDPMCSL_sieveCurve(minPar, maxPar, sieveCurveDiameter, sieveCurvePassing)
        else:
            newSieveCurveD, newSieveCurveP, w_min, w_max, NewSet = 0, 0, 0, 0, 0

        # Calculates volume of particles needed
        [volFracPar, parVolTotal, cdf, cdf1, kappa_i] = calc_LDPMCSL_parVolume(tetVolume, wcRatio, cementC,
                                                    airFrac, fullerCoef, 
                                                    flyashC, silicaC, scmC, fillerC,
                                                    flyashDensity, silicaDensity, 
                                                    scmDensity, fillerDensity, cementDensity,
                                                    densityWater, minPar, maxPar,
                                                    newSieveCurveD, newSieveCurveP, 
                                                    NewSet, w_min, w_max)



        print("Status: Calculating list of particles.") 
        # Calculate list of particle diameters for placement
        [maxParNum,parDiameterList] = gen_LDPMCSL_particleList(parVolTotal,minPar,maxPar,newSieveCurveD,\
            cdf,kappa_i,NewSet,fullerCoef)

        # Initialize empty particle nodes list outside geometry
        internalNodes = (np.zeros((len(parDiameterList),3))+2)*maxC













    ########################## Begin Placing Particles ##############################




    print('Status: Placing particles into geometry. (' + str(0) + '/' + str(len(internalNodes)) + ')') 

    # Initialize values
    newMaxIter = 6
    particlesPlaced = 0




    if multiMatToggle == "On":

        for i in range(3):

            # Place in order of aggregate > ITZ > binder
            if i == 0:
                [grainsDiameterList,voxels,grainMin,grainMax] = [aggGrainsDiameterList,aggVoxels,grainAggMin,grainAggMax]
            elif i == 1:
                [grainsDiameterList,voxels,grainMin,grainMax] = [itzGrainsDiameterList,itzVoxels,grainITZMin,grainITZMax]
            elif i == 2:
                [grainsDiameterList,voxels,grainMin,grainMax] = [binderGrainsDiameterList,binderVoxels,grainBinderMin,grainBinderMax]


            # Generate particles for length of needed aggregate (not placed via MPI)
            for x in range(particlesPlaced,len(grainsDiameterList)):

                # Generate particle
                [newMaxIter,node,iterReq] = gen_LDPMCSL_subParticle(surfaceNodes,grainsDiameterList[x],meshVertices,meshTets,newMaxIter,maxIter,grainMin,grainMax,\
                    aggOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes,\
                    multiMatX,multiMatY,multiMatZ,multiMatRes,voxels,minC,maxC)


                if len(grainsDiameterList)<=1000:
                    # Update number particles placed every 1%
                    if x % np.rint(len(grainsDiameterList)/100) == 0:
                        print("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')
                elif len(grainsDiameterList)<=10000:
                    # Update number particles placed every 0.1%
                    if x % np.rint(len(grainsDiameterList)/1000) == 0:
                        print("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')
                else:
                    # Update number particles placed every 0.01%
                    if x % np.rint(len(grainsDiameterList)/10000) == 0:
                        print("Status: Placing material " + str(i) + " grains into geometry. (" + str(x) + '/' + str(len(grainsDiameterList)) + ')')

                if i == 0:
                    internalNodes[x,:] = node
                elif i == 1:
                    internalNodes[x+len(aggGrainsDiameterList),:] = node
                elif i == 2:
                    internalNodes[x+len(aggGrainsDiameterList)+len(itzGrainsDiameterList),:] = node



            print("Status: Placing material " + str(i) + " grains into geometry. (" + str(len(grainsDiameterList)) + '/' + str(len(grainsDiameterList)) + ')')

        materialList = np.concatenate((np.ones(len(aggGrainsDiameterList))*3,np.ones(len(itzGrainsDiameterList))*1, np.ones(len(binderGrainsDiameterList))*2))

        # Set minimum particle to be smallest of the three materials 
        minPar = min(grainAggMin,grainITZMin,grainBinderMin)

    # Create empty lists if not cementStructure
    PoresDiameterList, ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0






    if multiMatToggle == "Off":

        if numCPU > 1:
        
            
            for increment in range(numIncrements-1):

                process_pool = multiprocessing.Pool(numCPU)

                outputMPI = process_pool.map(functools.partial(gen_LDPMCSL_particleMPI, surfaceNodes,maxParNum, minC, maxC, meshVertices, \
                    meshTets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minPar,\
                    maxPar,aggOffset,verbose,parDiameterList,maxEdgeLength,max_dist,internalNodes), parDiameterList[particlesPlaced:particlesPlaced+math.floor(len(parDiameterList)/numIncrements)])

                nodeMPI = np.array(outputMPI)[:,0:3]
                diameter = np.array(outputMPI)[:,3]
                newMaxIter = int(max(np.array(outputMPI)[:,4]))
                maxAttempts = int(max(np.array(outputMPI)[:,5]))

                particlesPlaced = particlesPlaced+len(np.array(outputMPI)[:,0:3])        

                for x in range(len(nodeMPI)):

                    # Store placed particles from this increment
                    internalNodes[particlesPlaced+x,:] = nodeMPI[x,:]

                    # Obtain extents for floating bin for node to test
                    binMin = np.array(([nodeMPI[x,0]-diameter[x]/2-maxPar/2-aggOffset,\
                        nodeMPI[x,1]-diameter[x]/2-maxPar/2-aggOffset,nodeMPI[x,2]-\
                        diameter[x]/2-maxPar/2-aggOffset]))
                    binMax = np.array(([nodeMPI[x,0]+diameter[x]/2+maxPar/2+aggOffset,\
                        nodeMPI[x,1]+diameter[x]/2+maxPar/2+aggOffset,nodeMPI[x,2]+\
                        diameter[x]/2+maxPar/2+aggOffset]))

                    # Check if particle overlapping any just added particles (ignore first one placed)
                    if x > 0:

                        overlap = check_LDPMCSL_particleOverlapMPI(nodeMPI[x,:],diameter[x],binMin,\
                            binMax,minPar,aggOffset,nodeMPI[0:x],diameter[0:x])

                        if overlap == True:

                            [newMaxIter,node,iterReq] = gen_LDPMCSL_particle(surfaceNodes,\
                                parDiameterList[particlesPlaced+x], meshVertices, \
                                meshTets,newMaxIter,maxIter,minPar,\
                                maxPar,aggOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes)
                            
                            internalNodes[particlesPlaced+x,:] = node[0,:]


                print("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')


        # Generate particles for length of needed aggregate (not placed via MPI)
        for x in range(particlesPlaced,len(parDiameterList)):

            # Generate particle
            [newMaxIter,node,iterReq] = gen_LDPMCSL_particle(surfaceNodes,parDiameterList[x],meshVertices,meshTets,newMaxIter,maxIter,minPar,maxPar,\
                aggOffset,parDiameterList,coord1,coord2,coord3,coord4,maxEdgeLength,max_dist,internalNodes)



            if len(parDiameterList)<=1000:
                # Update number particles placed every 1%
                if x % np.rint(len(parDiameterList)/100) == 0:
                    print("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')
            elif len(parDiameterList)<=10000:
                # Update number particles placed every 0.1%
                if x % np.rint(len(parDiameterList)/1000) == 0:
                    print("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')
            else:
                # Update number particles placed every 0.01%
                if x % np.rint(len(parDiameterList)/10000) == 0:
                    print("Status: Placing particles into geometry. (" + str(x) + '/' + str(len(parDiameterList)) + ')')

            internalNodes[x,:] = node

        print("Status: Placing particles into geometry. (" + str(len(parDiameterList)) + '/' + str(len(parDiameterList)) + ')')


        materialList = np.ones(len(parDiameterList))

        # Create empty lists if not multi-material or cementStructure
        aggGrainsDiameterList, itzGrainsDiameterList, binderGrainsDiameterList, PoresDiameterList,\
            ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0,0,0,0


    # Save the internalNodes list to a temporary file
    np.save(tempPath + 'internalNodes.npy', internalNodes)

    # Save the materialList to a temporary file
    np.save(tempPath + 'materialList.npy', materialList)

    # Save the parDiameterList to a temporary file
    np.save(tempPath + 'parDiameterList.npy', parDiameterList)

    if multiMatToggle == "Off":

        # Save the volFracPar to a temporary file
        np.save(tempPath + 'volFracPar.npy', volFracPar)




    #placementTime = round(time.time() - start_time,2)   
    nParticles = len(parDiameterList)