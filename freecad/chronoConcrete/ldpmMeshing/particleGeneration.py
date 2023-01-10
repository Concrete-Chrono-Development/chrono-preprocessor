#!usr/bin/python
# -*- coding: utf-8 -*-

# Import required packages
import sys
import os
import re
import shutil
import math
import numpy as np
import time
import subprocess
from pathlib import Path
import pickle
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as mticker
from scipy import spatial
from array import array
import scipy
from scipy.integrate import quad
from numpy import linspace
import multiprocessing
import functools

# Initialize code start time to measure performance
start_time = time.time()

# Increase system check interval to improve performance
sys.setcheckinterval(1000)

# Turn off error for divide by zero and invalid operations
np.seterr(divide='ignore', invalid='ignore')


####################################################################################################
## Main Class Operations                                                                          ##
####################################################################################################


class ParticleGen:
    def __init__(self, maxAggD,minAggD,fullerCoef,wcRatio,cementC,\
        volFracAir,q,maxIter,geoFile,aggOffsetCoeff,densityWater,densityCement,\
        dataType,output,verbose,fibers,dFiber,lFiber,vFiber,fiberFile,\
        multiMaterial,materialFile,maxGrainD,minGrainD,grainFullerCoef,\
        maxBinderD,minBinderD,binderFullerCoef,maxITZD,minITZD,ITZFullerCoef,\
        rebar,rebarFile1,dRebar1,rebarFile2,dRebar2,rebarFile3,dRebar3,edgeElements,\
        surfaceExtLength,fiberOrientation,orientationStrength,sieveCurveDiameter,sieveCurvePassing,\
        grainSieveCurveDiameter,grainSieveCurvePassing,binderSieveCurveDiameter,\
        binderSieveCurvePassing,itzSieveCurveDiameter,itzSieveCurvePassing,materialRule,\
        cementStructure,cementmaterialFile,cementStructureResolution,numberofPhases,\
        maxPoresD,minPoresD,PoresFullerCoef,PoresSieveCurveDiameter,PoresSieveCurvePassing,\
        maxClinkerD,minClinkerD,ClinkerFullerCoef,ClinkerSieveCurveDiameter,ClinkerSieveCurvePassing,\
        maxCHD,minCHD,CHFullerCoef,CHSieveCurveDiameter,CHSieveCurvePassing,\
        maxCSH_LDD,minCSH_LDD,CSH_LDFullerCoef,CSH_LDSieveCurveDiameter,CSH_LDSieveCurvePassing,\
        maxCSH_HDD,minCSH_HDD,CSH_HDFullerCoef,CSH_HDSieveCurveDiameter,CSH_HDSieveCurvePassing,\
        Mean_CSH_LD,Mean_CSH_HD,SDev_CSH_LD,SDev_CSH_HD,Alphac,SaturatedCSHDensity,\
        caeFile,periodic,prismX,prismY,prismZ,randomField,cutFiber,outputUnits,numIncrements,numCPU):
        
        # Convert density to Kg/m3

        cementC = cementC*1.0E+12
        densityCement = densityCement*1.0E+12
        densityWater = densityWater*1.0E+12


        if periodic in ['on','On','Y','y','Yes','yes']:
            geoFile = 'prism_' + str(prismX).replace('.', '') + 'mmx' + str(prismY).replace('.', '') + 'mmx' + \
                str(prismZ).replace('.', '') + 'mm.mesh'

        # Basic Calcs
        aggOffset = aggOffsetCoeff*minAggD

        
        if caeFile :
            cur_path = os.path.join(os.getcwd(), 'geometry')
            cmd = 'abaqus cae noGUI=cae2vtk.py -- "{}" "{}"'.format(caeFile, cur_path)
            print(cmd)
            subprocess.check_output(cmd, shell=True)
            caeFile = os.path.basename(caeFile)
            caeName, caeExtension = caeFile.split('.')
            geoFile = caeName + '.vtk'
            for item in os.listdir(os.getcwd()):
                if "abaqus.rpy" in item or "abaqus_acis.log" in item:
                    os.remove(os.path.join(os.getcwd(), item))

        geoName = geoFile.split(".")[0]
        geoExtension = geoFile.split(".")[1]

        # Remove directory/files if exists already
        try:
            shutil.rmtree(Path('meshes/' + geoName))
        except:
            pass

        # try:
        #     shutil.rmtree(Path('temp'))
        # except:
        #     pass

        if not os.path.exists(Path('temp')):
            os.mkdir(Path('temp'))

        if not os.path.exists(Path('meshes/' + geoName)):
            os.mkdir(Path('meshes/' + geoName))

        if periodic in ['on','On','Y','y','Yes','yes']:
            periodicMesh = self.periodicMesher(prismX,prismY,prismZ,minAggD,\
                maxAggD,geoName,geoFile)


        # Else mesh the concrete geometry
        else:
            concreteMesh = self.meshConcreteGeometry(geoName,minAggD,2*minAggD,\
                geoFile,geoExtension)


        if rebar in ['on','On','Y','y','Yes','yes']:
            # Mesh the rebar
            rebarData = self.meshRebarGeometry(minAggD,maxAggD,rebar,rebarFile1,\
                dRebar1,rebarFile2,dRebar2,rebarFile3,dRebar3,geoName)
        else:
            pass

        print('Starting particle placement.')

        # Reformat the mesh into array of nodes and tets
        command1 = Path('temp/' + geoName + '.mesh')
        command2 = Path('temp/' + geoName + '2D.mesh')
        vertices = self.readConcreteMesh(command1,command2)

        # Gets extents of geometry
        [minC,maxC] = self.meshExtents(vertices)

        # Gets volume of geometry
        tetVolume = self.meshVolume(vertices,self.tets)

        # Generates points for all external triangles
        trianglePoints = self.formTriangles(vertices, self.triangles)

        # Initialize list of failed placement trials (improve speed by 
        # immediately discarding trials which overlap these points)
        #ParticleGen.badItem = 0
        #ParticleGen.badList = np.array([maxC,]*100)*2

        tets = (self.tets).astype(int)
        
        coord1 = vertices[tets[:,0]-1]
        coord2 = vertices[tets[:,1]-1]
        coord3 = vertices[tets[:,2]-1]
        coord4 = vertices[tets[:,3]-1]

        if multiMaterial in ['on','On','Y','y','Yes','yes']:
            
            # Read the voxelated microstructure file and store
            [X_Size, Y_Size, Z_Size, resolution, voxels] = self.readMultiMaterial(materialFile)

            # Check if the voxelated microstructure is larger than the geometry
            sizeCheck = self.checkMicrostructureSize(X_Size, Y_Size, Z_Size, \
                resolution, minC, maxC)

            # Organize and store voxels of each material
            [agg,itz,binder] = self.sortMicrostructure(voxels)

            # Aggregate Grain Sieve Curve
            if grainSieveCurveDiameter != []:
                # Shifts sieve curve to appropriate range
                [newGrainSieveCurveD,newGrainSieveCurveP,grainNewSet,grainW_min,grainW_max] = \
                    self.sieveCurve(minGrainD,maxGrainD,grainSieveCurveDiameter,grainSieveCurvePassing)
            else:
                newGrainSieveCurveD, newGrainSieveCurveP, grainW_min, grainW_max, grainNewSet = 0, 0, 0, 0, 0

            # Binder Grain Sieve Curve
            if binderSieveCurveDiameter != []:
                # Shifts sieve curve to appropriate range
                [newBinderSieveCurveD,newBinderSieveCurveP,binderNewSet,binderW_min,binderW_max] = \
                    self.sieveCurve(minBinderD,maxBinderD,binderSieveCurveDiameter,binderSieveCurvePassing)
            else:
                newBinderSieveCurveD, newBinderSieveCurveP, binderW_min, binderW_max, binderNewSet = 0, 0, 0, 0, 0

            # ITZ Grain Sieve Curve
            if itzSieveCurveDiameter != []:
                # Shifts sieve curve to appropriate range
                [newITZSieveCurveD,newITZSieveCurveP,itzNewSet,itzW_min,itzW_max] = self.sieveCurve(minITZD,maxITZD,\
                    itzSieveCurveDiameter,itzSieveCurvePassing)
            else:
                newITZSieveCurveD, newITZSieveCurveP, itzW_min, itzW_max, itzNewSet = 0, 0, 0, 0, 0

            # Calculates volume of aggregate grains needed
            [volAggGrains,cdfAgg,cdf1Agg,kappa_iAgg] = self.aggVolume(tetVolume*len(agg)/(len(agg)+\
                len(itz)+len(binder)),wcRatio,cementC,volFracAir,\
                grainFullerCoef,densityCement,densityWater,minGrainD,maxGrainD,\
                newGrainSieveCurveD,newGrainSieveCurveP,grainNewSet,grainW_min,grainW_max)

            # Calculates volume of binder particles needed
            [volBinder,cdfBinder,cdf1Binder,kappa_iBinder] = self.aggVolume(tetVolume*len(binder)/(len(agg)+\
                len(itz)+len(binder)),wcRatio,cementC,\
                volFracAir,binderFullerCoef,densityCement,densityWater,\
                minBinderD,maxBinderD,newBinderSieveCurveD,newBinderSieveCurveP,binderNewSet,binderW_min,binderW_max)

            # Calculates volume of ITZ particles needed
            [volITZ,cdfITZ,cdf1ITZ,kappa_iITZ] = self.aggVolume(tetVolume*len(itz)/(len(agg)+len(itz)+\
                len(binder)),wcRatio,cementC,\
                volFracAir,ITZFullerCoef,densityCement,densityWater,minITZD,maxITZD,\
                newITZSieveCurveD,newITZSieveCurveP,itzNewSet,itzW_min,itzW_max)

            # Generates list of needed aggregate grains
            [maxAggGrainsNum,aggGrainsDiameterList] = self.aggList(\
                volAggGrains,minGrainD,maxGrainD,3-grainFullerCoef,newGrainSieveCurveD,cdfAgg,kappa_iAgg,grainNewSet)

            # Generates list of needed binder particles
            [maxBinderNum,binderDiameterList] = self.aggList(volBinder,\
                minBinderD,maxBinderD,3-binderFullerCoef,newBinderSieveCurveD,cdfBinder,kappa_iBinder,binderNewSet)

            # Generates list of needed ITZ particles
            [maxITZNum,itzDiameterList] = self.aggList(volITZ,minITZD,maxITZD,3-ITZFullerCoef,newITZSieveCurveD,cdfITZ,kappa_iITZ,itzNewSet)

            # Initialize empty aggregate grain nodes list outside geometry
            self.nodesAgg = (np.zeros((len(aggGrainsDiameterList),3))+2)*maxC

            # Initialize empty binder particle nodes list outside geometry
            self.nodesBinder = (np.zeros((len(binderDiameterList),3))+2)*maxC

            # Initialize empty ITZ particle nodes list outside geometry
            self.nodesITZ = (np.zeros((len(itzDiameterList),3))+2)*maxC
            
            # Initialize empty list of all nodes outside geometry
            self.nodes = (np.zeros((len(aggGrainsDiameterList)+\
                len(binderDiameterList)+len(itzDiameterList),3))+2)*maxC

            # Calculate surface mesh max edge length
            maxEdgeLength1 = max(np.linalg.norm(self.vertices2D[self.triangles2D[:,1].astype(int)-1,0:3]-self.vertices2D[self.triangles2D[:,0].astype(int)-1,0:3], axis=1))
            maxEdgeLength2 = max(np.linalg.norm(self.vertices2D[self.triangles2D[:,2].astype(int)-1,0:3]-self.vertices2D[self.triangles2D[:,1].astype(int)-1,0:3], axis=1))
            maxEdgeLength3 = max(np.linalg.norm(self.vertices2D[self.triangles2D[:,2].astype(int)-1,0:3]-self.vertices2D[self.triangles2D[:,0].astype(int)-1,0:3], axis=1))
            maxEdgeLength = max([maxEdgeLength1,maxEdgeLength2,maxEdgeLength3])


            print('Aggregate grains...')

            # Generate grains for length of needed aggregate grains
            for x in range(0,len(aggGrainsDiameterList)):

                # Generate particle
                genParticles = self.generateSubParticle(x,trianglePoints,\
                    aggGrainsDiameterList[x],maxAggGrainsNum, minC, maxC, \
                    vertices,self.tets, coord1,coord2,coord3,coord4,maxIter,\
                    minGrainD,maxGrainD,aggOffsetCoeff*minGrainD,verbose,\
                    aggGrainsDiameterList,agg,X_Size,Y_Size,Z_Size,resolution,\
                    aggGrainsDiameterList, binderDiameterList, itzDiameterList,\
                    'Aggregate',maxEdgeLength)
                
                self.nodesAgg[x,:] = self.node[0,:]
                self.nodes[x,:] = self.node[0,:]

            print('ITZ particles...')

            # Generate particles for length of needed binder particles
            for x in range(0,len(itzDiameterList)):

                # Generate particle
                genParticles = self.generateSubParticle(x,trianglePoints,\
                    itzDiameterList[x],maxITZNum, minC, maxC, vertices, \
                    self.tets, coord1,coord2,coord3,coord4,maxIter,minITZD,\
                    maxITZD,aggOffsetCoeff*minITZD,verbose,itzDiameterList,itz,X_Size,\
                    Y_Size,Z_Size,resolution,aggGrainsDiameterList,\
                    binderDiameterList, itzDiameterList, 'ITZ',maxEdgeLength)
                
                self.nodesITZ[x,:] = self.node[0,:]
                self.nodes[x+len(aggGrainsDiameterList),:] = self.node[0,:]

            print('Binder particles...')

            # Generate particles for length of needed binder particles
            for x in range(0,len(binderDiameterList)):

                # Generate particle
                genParticles = self.generateSubParticle(x,trianglePoints,\
                    binderDiameterList[x],maxBinderNum, minC, maxC, vertices, \
                    self.tets, coord1,coord2,coord3,coord4,maxIter,minBinderD,maxBinderD,\
                    aggOffsetCoeff*minBinderD,verbose,binderDiameterList,binder,X_Size,Y_Size,Z_Size,resolution,\
                    aggGrainsDiameterList, binderDiameterList, itzDiameterList, 'Binder',maxEdgeLength)
                
                self.nodesBinder[x,:] = self.node[0,:]
                self.nodes[x+len(aggGrainsDiameterList)+len(itzDiameterList),:] = self.node[0,:]

            aggDiameterList = np.concatenate((aggGrainsDiameterList, itzDiameterList, binderDiameterList))
            materialList = np.concatenate((np.ones(len(aggGrainsDiameterList))*3,\
                np.ones(len(itzDiameterList))*1, np.ones(len(binderDiameterList))*2))
            minAggD = minBinderD

            placementTime = round(time.time() - start_time,2)   
            nParticles = len(aggGrainsDiameterList)+len(binderDiameterList)\
                +len(itzDiameterList)

            print(str(nParticles) + ' particles/grains placed in ' + \
                str(placementTime) + ' seconds' )
            print('----------------------------------------------------')

            # Create empty lists if not cementStructure
            PoresDiameterList, ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0

            print(str(nParticles) + ' particles/grains placed in ' + \
                str(placementTime) + ' seconds' )
            print('----------------------------------------------------')

        elif cementStructure in ['on','On','Y','y','Yes','yes']:

            # Read the voxelated cementStructure file and store
            [CementStructure_X_Size, CementStructure_Y_Size, CementStructure_Z_Size,\
                XYZID, CementStructureVoxels] = self.readCementStructure(cementmaterialFile,minC,cementStructureResolution)

            # Check if the voxelated cementStructure is larger than the geometry
            sizeCheck = self.checkCementStructureSize(CementStructure_X_Size,CementStructure_Y_Size,\
                CementStructure_Z_Size,cementStructureResolution,minC,maxC)

            # Organize and store voxels of each material
            [PoresVoxelData, ClinkerVoxelData, CHVoxelData, CSH_LDVoxelData,\
                CSH_HDVoxelData, AllHydrationProduct] = self.sortCementStructure(CementStructureVoxels,\
                    XYZID,wcRatio,Alphac,SaturatedCSHDensity,Mean_CSH_HD,Mean_CSH_LD,SDev_CSH_HD,SDev_CSH_LD,densityWater)

            # Capillary Pores Sieve Curve
            if PoresSieveCurveDiameter != []:
                # Shifts sieve curve to appropriate range
                [newPoresSieveCurveD,newPoresSieveCurveP,PoresNewSet,PoresW_min,PoresW_max] = \
                    self.sieveCurve(minPoresD,maxPoresD,PoresSieveCurveDiameter,PoresSieveCurvePassing)
            else:
                newPoresSieveCurveD, newPoresSieveCurveP, PoresNewSet, PoresW_min, PoresW_max = 0, 0, 0, 0, 0

            # Clinker Sieve Curve
            if ClinkerSieveCurveDiameter != []:
                # Shifts sieve curve to appropriate range
                [newClinkerSieveCurveD,newClinkerSieveCurveP,ClinkerNewSet,ClinkerW_min,ClinkerW_max] = \
                    self.sieveCurve(minClinkerD,maxClinkerD,ClinkerSieveCurveDiameter,ClinkerSieveCurvePassing)
            else:
                newClinkerSieveCurveD, newClinkerSieveCurveP, ClinkerNewSet, ClinkerW_min, ClinkerW_max = 0, 0, 0, 0, 0

            # CH Sieve Curve
            if CHSieveCurveDiameter != []:
                # Shifts sieve curve to appropriate range
                [newCHSieveCurveD,newCHSieveCurveP,CHNewSet,CHW_min,CHW_max] = \
                    self.sieveCurve(minCHD,maxCHD,CHSieveCurveDiameter,CHSieveCurvePassing)
            else:
                newCHSieveCurveD, newCHSieveCurveP, CHNewSet, CHW_min, CHW_max = 0, 0, 0, 0, 0

            # CSH_LD Sieve Curve
            if CSH_LDSieveCurveDiameter != []:
                # Shifts sieve curve to appropriate range
                [newCSH_LDSieveCurveD,newCSH_LDSieveCurveP,CSH_LDNewSet,CSH_LDW_min,CSH_LDW_max] = \
                    self.sieveCurve(minCSH_LDD,maxCSH_LDD,CSH_LDSieveCurveDiameter,CSH_LDSieveCurvePassing)
            else:
                newCSH_LDSieveCurveD, newCSH_LDSieveCurveP, CSH_LDNewSet, CSH_LDW_min, CSH_LDW_max = 0, 0, 0, 0, 0  

            # CSH_HD Sieve Curve
            if CSH_HDSieveCurveDiameter != []:
                # Shifts sieve curve to appropriate range
                [newCSH_HDSieveCurveD,newCSH_HDSieveCurveP,CSH_HDNewSet,CSH_HDW_min,CSH_HDW_max] = \
                    self.sieveCurve(minCSH_HDD,maxCSH_HDD,CSH_HDSieveCurveDiameter,CSH_HDSieveCurvePassing)
            else:
                newCSH_HDSieveCurveD, newCSH_HDSieveCurveP, CSH_HDNewSet, CSH_HDW_min, CSH_HDW_max = 0, 0, 0, 0, 0

            # Calculates volume of Capillary Pores needed
            [volPores,cdfPores,cdf1Pores,kappa_iPores] = self.aggVolume(tetVolume*len(PoresVoxelData)/\
                len(AllHydrationProduct),wcRatio,cementC,volFracAir,PoresFullerCoef,densityCement,\
                densityWater,minPoresD,maxPoresD,newPoresSieveCurveD,newPoresSieveCurveP,PoresNewSet,PoresW_min,PoresW_max)

            # Calculates volume of Clinker needed
            [volClinker,cdfClinker,cdf1Clinker,kappa_iClinker] = self.aggVolume(tetVolume*len(ClinkerVoxelData)/\
                len(AllHydrationProduct),wcRatio,cementC,volFracAir,ClinkerFullerCoef,densityCement,\
                densityWater,minClinkerD,maxClinkerD,newClinkerSieveCurveD,newClinkerSieveCurveP,ClinkerNewSet,ClinkerW_min,ClinkerW_max)

            # Calculates volume of CH needed
            [volCH,cdfCH,cdf1CH,kappa_iCH] = self.aggVolume(tetVolume*len(CHVoxelData)/\
                len(AllHydrationProduct),wcRatio,cementC,volFracAir,CHFullerCoef,densityCement,\
                densityWater,minCHD,maxCHD,newCHSieveCurveD,newCHSieveCurveP,CHNewSet,CHW_min,CHW_max)

            # Calculates volume of CSH_LD needed
            [volCSH_LD,cdfCSH_LD,cdf1CSH_LD,kappa_iCSH_LD] = self.aggVolume(tetVolume*len(CSH_LDVoxelData)/\
                len(AllHydrationProduct),wcRatio,cementC,volFracAir,CSH_LDFullerCoef,densityCement,\
                densityWater,minCSH_LDD,maxCSH_LDD,newCSH_LDSieveCurveD,newCSH_LDSieveCurveP,CSH_LDNewSet,CSH_LDW_min,CSH_LDW_max)

            # Calculates volume of CSH_HD needed
            [volCSH_HD,cdfCSH_HD,cdf1CSH_HD,kappa_iCSH_HD] = self.aggVolume(tetVolume*len(CSH_HDVoxelData)/\
                len(AllHydrationProduct),wcRatio,cementC,volFracAir,CSH_HDFullerCoef,densityCement,\
                densityWater,minCSH_HDD,maxCSH_HDD,newCSH_HDSieveCurveD,newCSH_HDSieveCurveP,CSH_HDNewSet,CSH_HDW_min,CSH_HDW_max)

            # Generates list of needed Capillary Pores
            [maxPoresNum,PoresDiameterList] = self.aggList(volPores,\
                minPoresD,maxPoresD,3-PoresFullerCoef,newPoresSieveCurveD,cdfPores,kappa_iPores,PoresNewSet)

            # Generates list of needed Clinker
            [maxClinkerNum,ClinkerDiameterList] = self.aggList(volClinker,\
                minClinkerD,maxClinkerD,3-ClinkerFullerCoef,newClinkerSieveCurveD,cdfClinker,kappa_iClinker,ClinkerNewSet)

            # Generates list of needed CH
            [maxCHNum,CHDiameterList] = self.aggList(volCH,\
                minCHD,maxCHD,3-CHFullerCoef,newCHSieveCurveD,cdfCH,kappa_iCH,CHNewSet)

            # Generates list of needed CSH_LD
            [maxCSH_LDNum,CSH_LDDiameterList] = self.aggList(volCSH_LD,\
                minCSH_LDD,maxCSH_LDD,3-CSH_LDFullerCoef,newCSH_LDSieveCurveD,cdfCSH_LD,kappa_iCSH_LD,CSH_LDNewSet)

            # Generates list of needed CSH_HD
            [maxCSH_HDNum,CSH_HDDiameterList] = self.aggList(volCSH_HD,\
                minCSH_HDD,maxCSH_HDD,3-CSH_HDFullerCoef,newCSH_HDSieveCurveD,cdfCSH_HD,kappa_iCSH_HD,CSH_HDNewSet)

            # Initialize empty Capillary Pores particle nodes list outside geometry
            self.nodesPores = (np.zeros((len(PoresDiameterList),3))+2)*maxC

            # Initialize empty Clinker particle nodes list outside geometry
            self.nodesClinker = (np.zeros((len(ClinkerDiameterList),3))+2)*maxC

            # Initialize empty CH particle nodes list outside geometry
            self.nodesCH = (np.zeros((len(CHDiameterList),3))+2)*maxC

            # Initialize empty CSH_LD particle nodes list outside geometry
            self.nodesCSH_LD = (np.zeros((len(CSH_LDDiameterList),3))+2)*maxC

            # Initialize empty CSH_HD particle nodes list outside geometry
            self.nodesCSH_HD = (np.zeros((len(CSH_HDDiameterList),3))+2)*maxC

            # Initialize empty list of all nodes outside geometry
            self.nodes = (np.zeros((len(PoresDiameterList)+\
                len(ClinkerDiameterList)+len(CHDiameterList)+
                len(CSH_LDDiameterList)+len(CSH_HDDiameterList),3))+2)*maxC

            print('Capillary Pores...')

            # Generate grains for length of needed Capillary Pores
            for x in range(0,len(PoresDiameterList)):

                # Generate particle
                genParticles = self.generateSubParticleCement(x,trianglePoints,\
                    PoresDiameterList[x],maxPoresNum, minC, maxC, vertices, \
                    self.tets, coord1,coord2,coord3,coord4,maxIter,minPoresD,\
                    maxPoresD,aggOffsetCoeff*minPoresD,verbose,PoresDiameterList,\
                    PoresVoxelData[:,0],CementStructure_X_Size,CementStructure_Y_Size,\
                    CementStructure_Z_Size,cementStructureResolution,PoresDiameterList,\
                    ClinkerDiameterList,CHDiameterList,CSH_LDDiameterList,CSH_HDDiameterList,'Pores',maxEdgeLength)
                
                self.nodesPores[x,:] = self.node[0,:]
                self.nodes[x,:] = self.node[0,:]

            print('Clinker...')

            # Generate grains for length of needed Clinker
            for x in range(0,len(ClinkerDiameterList)):

                # Generate particle
                genParticles = self.generateSubParticleCement(x,trianglePoints,\
                    ClinkerDiameterList[x],maxClinkerNum, minC, maxC, vertices, \
                    self.tets, coord1,coord2,coord3,coord4,maxIter,minClinkerD,\
                    maxClinkerD,aggOffsetCoeff*minClinkerD,verbose,ClinkerDiameterList,\
                    ClinkerVoxelData[:,0],CementStructure_X_Size,CementStructure_Y_Size,\
                    CementStructure_Z_Size,cementStructureResolution,PoresDiameterList,\
                    ClinkerDiameterList,CHDiameterList,CSH_LDDiameterList,CSH_HDDiameterList, 'Clinker',maxEdgeLength)
                
                self.nodesClinker[x,:] = self.node[0,:]
                self.nodes[x+len(PoresDiameterList),:] = self.node[0,:]

            # Generate grains for length of needed CH
            for x in range(0,len(CHDiameterList)):

                # Generate particle
                genParticles = self.generateSubParticleCement(x,trianglePoints,\
                    CHDiameterList[x],maxCHNum, minC, maxC, vertices, \
                    self.tets, coord1,coord2,coord3,coord4,maxIter,minCHD,\
                    maxCHD,aggOffsetCoeff*minCHD,verbose,CHDiameterList,\
                    CHVoxelData[:,0],CementStructure_X_Size,CementStructure_Y_Size,\
                    CementStructure_Z_Size,cementStructureResolution,PoresDiameterList,\
                    ClinkerDiameterList,CHDiameterList,CSH_LDDiameterList,CSH_HDDiameterList, 'CH',maxEdgeLength)
                
                self.nodesCH[x,:] = self.node[0,:]
                self.nodes[x+len(PoresDiameterList)+len(ClinkerDiameterList),:] = self.node[0,:]

            # Generate grains for length of needed CSH_LD
            for x in range(0,len(CSH_LDDiameterList)):

                # Generate particle
                genParticles = self.generateSubParticleCement(x,trianglePoints,\
                    CSH_LDDiameterList[x],maxCSH_LDNum, minC, maxC, vertices, \
                    self.tets, coord1,coord2,coord3,coord4,maxIter,minCSH_LDD,\
                    maxCSH_LDD,aggOffsetCoeff*minCSH_LDD,verbose,CSH_LDDiameterList,\
                    CSH_LDVoxelData[:,0],CementStructure_X_Size,CementStructure_Y_Size,\
                    CementStructure_Z_Size,cementStructureResolution,PoresDiameterList,\
                    ClinkerDiameterList,CHDiameterList,CSH_LDDiameterList,CSH_HDDiameterList, 'CSH_LD',maxEdgeLength)
                
                self.nodesCSH_LD[x,:] = self.node[0,:]
                self.nodes[x+len(PoresDiameterList)+len(ClinkerDiameterList)+len(CHDiameterList),:] = self.node[0,:]

            # Generate grains for length of needed CSH_HD
            for x in range(0,len(CSH_HDDiameterList)):

                # Generate particle
                genParticles = self.generateSubParticleCement(x,trianglePoints,\
                    CSH_HDDiameterList[x],maxCSH_HDNum, minC, maxC, vertices, \
                    self.tets, coord1,coord2,coord3,coord4,maxIter,minCSH_HDD,\
                    maxCSH_HDD,aggOffsetCoeff*minCSH_HDD,verbose,CSH_HDDiameterList,\
                    CSH_HDVoxelData[:,0],CementStructure_X_Size,CementStructure_Y_Size,\
                    CementStructure_Z_Size,cementStructureResolution,PoresDiameterList,\
                    ClinkerDiameterList,CHDiameterList,CSH_LDDiameterList,CSH_HDDiameterList, 'CSH_HD',maxEdgeLength)
                
                self.nodesCSH_HD[x,:] = self.node[0,:]
                self.nodes[x+len(PoresDiameterList)+len(ClinkerDiameterList)+len(CHDiameterList)+len(CSH_LDDiameterList),:] = self.node[0,:]

            aggDiameterList = np.concatenate((PoresDiameterList, ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList))
            materialList = np.concatenate((np.ones(len(PoresDiameterList))*1,\
                np.ones(len(ClinkerDiameterList))*2, np.ones(len(CHDiameterList))*3,\
                np.ones(len(CSH_LDDiameterList))*4, np.ones(len(CSH_HDDiameterList))*5))
            print(materialList)
            print(materialList.shape)
            minAggD = min(minPoresD,minClinkerD,minCHD,minCSH_LDD,minCSH_HDD)
            print(minAggD)

            # Generates list of needed CSH_HD
            [maxCSH_HDNum,CSH_HDDiameterList] = self.aggList(volCSH_HD,\
                minCSH_HDD,maxCSH_HDD,3-CSH_HDFullerCoef,newCSH_HDSieveCurveD,cdfCSH_HD,kappa_iCSH_HD,CSH_HDNewSet)

            placementTime = round(time.time() - start_time,2)   
            nParticles = len(PoresDiameterList)+len(ClinkerDiameterList)\
                +len(CHDiameterList)+len(CSH_LDDiameterList)+\
                len(CSH_HDDiameterList)

            # Create empty lists if not multi-material
            aggGrainsDiameterList, itzDiameterList, binderDiameterList = 0,0,0

            print(str(nParticles) + ' particles/ placed in ' + \
                str(placementTime) + ' seconds' )
            print('----------------------------------------------------')      

        else:

            if sieveCurveDiameter != []:

                # Shifts sieve curve to appropriate range
                [newSieveCurveD,newSieveCurveP,NewSet,w_min,w_max] = self.sieveCurve(minAggD,maxAggD,\
                    sieveCurveDiameter,sieveCurvePassing)

            else:

                newSieveCurveD, newSieveCurveP, w_min, w_max, NewSet = 0, 0, 0, 0, 0

            # Calculates volume of aggregate needed
            [aggVolTotal,cdf,cdf1,kappa_i] = self.aggVolume(tetVolume,wcRatio,cementC,\
                volFracAir,fullerCoef,densityCement,densityWater,minAggD,maxAggD,\
                newSieveCurveD,newSieveCurveP,NewSet,w_min,w_max)

            # Generates list of needed aggreate
            [maxAggNum,aggDiameterList] = self.aggList(aggVolTotal,minAggD,maxAggD,q,newSieveCurveD,\
                cdf,kappa_i,NewSet)

            # Initialize empty particle nodes list outside geometry
            self.nodes = (np.zeros((len(aggDiameterList),3))+2)*maxC

            # Calculate surface mesh max edge length
            maxEdgeLength1 = max(np.linalg.norm(self.vertices2D[self.triangles2D[:,1].astype(int)-1,0:3]-self.vertices2D[self.triangles2D[:,0].astype(int)-1,0:3], axis=1))
            maxEdgeLength2 = max(np.linalg.norm(self.vertices2D[self.triangles2D[:,2].astype(int)-1,0:3]-self.vertices2D[self.triangles2D[:,1].astype(int)-1,0:3], axis=1))
            maxEdgeLength3 = max(np.linalg.norm(self.vertices2D[self.triangles2D[:,2].astype(int)-1,0:3]-self.vertices2D[self.triangles2D[:,0].astype(int)-1,0:3], axis=1))
            maxEdgeLength = max([maxEdgeLength1,maxEdgeLength2,maxEdgeLength3])

            # Initialize values
            self.newMaxIter = 2
            newMaxIter = 2
            particlesPlaced = 0

            if numCPU > 1:
            
                if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
                    print("%s Remaining." % (len(aggDiameterList)))

                for increment in range(numIncrements-1):

                    process_pool = multiprocessing.Pool(numCPU)

                    outputMPI = process_pool.map(functools.partial(self.generateParticleMPI, trianglePoints,maxAggNum, minC, maxC, vertices, \
                        self.tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minAggD,\
                        maxAggD,aggOffset,verbose,aggDiameterList,maxEdgeLength), aggDiameterList[particlesPlaced:particlesPlaced+math.floor(len(aggDiameterList)/numIncrements)])

                    nodeMPI = np.array(outputMPI)[:,0:3]
                    diameter = np.array(outputMPI)[:,3]
                    newMaxIter = int(max(np.array(outputMPI)[:,4]))
                    maxAttempts = int(max(np.array(outputMPI)[:,5]))

                    particlesPlaced = particlesPlaced+len(np.array(outputMPI)[:,0:3])        

                    for x in range(len(nodeMPI)):

                        # Store placed particles from this increment
                        self.nodes[particlesPlaced+x,:] = nodeMPI[x,:]

                        # Obtain extents for floating bin for node to test
                        binMin = np.array(([nodeMPI[x,0]-diameter[x]/2-maxAggD/2-aggOffset,\
                            nodeMPI[x,1]-diameter[x]/2-maxAggD/2-aggOffset,nodeMPI[x,2]-\
                            diameter[x]/2-maxAggD/2-aggOffset]))
                        binMax = np.array(([nodeMPI[x,0]+diameter[x]/2+maxAggD/2+aggOffset,\
                            nodeMPI[x,1]+diameter[x]/2+maxAggD/2+aggOffset,nodeMPI[x,2]+\
                            diameter[x]/2+maxAggD/2+aggOffset]))

                        # Check if particle overlapping any just added particles (ignore first one placed)
                        if x > 0:

                            overlap = self.overlapCheckMPI(nodeMPI[x,:],diameter[x],binMin,\
                                binMax,minAggD,aggOffset,nodeMPI[0:x],diameter[0:x])

                            if overlap == True:

                                newMaxIter = self.generateParticle(particlesPlaced+x,trianglePoints,\
                                    aggDiameterList[particlesPlaced+x],maxAggNum, minC, maxC, vertices, \
                                    self.tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minAggD,\
                                    maxAggD,aggOffset,'No',aggDiameterList,maxEdgeLength)
                                
                                self.nodes[particlesPlaced+x,:] = self.node[0,:]

                    if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
                        print("%s Remaining. Maximum attempts required in increment: %s" % \
                            (len(aggDiameterList)-particlesPlaced, maxAttempts))


            # Generate particles for length of needed aggregate (not placed via MPI)
            for x in range(particlesPlaced,len(aggDiameterList)):

                # Generate particle
                newMaxIter = self.generateParticle(x,trianglePoints,\
                    aggDiameterList[x],maxAggNum, minC, maxC, vertices, \
                    self.tets, coord1,coord2,coord3,coord4,newMaxIter,maxIter,minAggD,\
                    maxAggD,aggOffset,verbose,aggDiameterList,maxEdgeLength)
                
                self.nodes[x,:] = self.node[0,:]
            


            materialList = np.ones(len(aggDiameterList))

            placementTime = round(time.time() - start_time,2)   
            nParticles = len(aggDiameterList)

            # Create empty lists if not multi-material or cementStructure
            aggGrainsDiameterList, itzDiameterList, binderDiameterList, PoresDiameterList,\
                ClinkerDiameterList, CHDiameterList, CSH_LDDiameterList, CSH_HDDiameterList = 0,0,0,0,0,0,0,0

            print(str(nParticles) + ' particles placed in ' + \
                str(placementTime) + ' seconds' )
            print('----------------------------------------------------')

        if fibers in ['on','On','Y','y','Yes','yes']:
            
            fiberStartTime = time.time()

            # Use fiber data from CT scan if file exists
            try:

                ctData = Path('fibers/' + fiberFile)

                # CTScanfiber data arrangmentt
                CTScanFiberData = self.CTScanFibers(ctData)
                CTScanFiberData = np.array(CTScanFiberData).reshape(-1,10)
                self.p1Fibers = CTScanFiberData[:,0:3]
                self.p2Fibers = CTScanFiberData[:,3:6]
                self.orienFibers = CTScanFiberData[:,6:9]
                self.fiberLengths = CTScanFiberData[:,9:10]


            # Generate fibers if no CT data
            except:

                print('Starting fibers placement.')

                if orientationStrength<0 or orientationStrength>1:
                    print('Fiber orientation strength is out of range, use 0-1')
                    print('Now exitting...')
                    exit()

                # Calculate number of fibers needed 
                nFiber = int(round(4*tetVolume*vFiber/(math.pi*dFiber**2*lFiber)))

                # Initialize empty fiber nodes list outside geometry
                self.p1Fibers = (np.zeros((nFiber,3))+2)*maxC
                self.p2Fibers = (np.zeros((nFiber,3))+2)*maxC
                self.orienFibers = (np.zeros((nFiber,3))+2)*maxC
                self.fiberLengths = (np.zeros((nFiber,1)))

                # Generate fibers for number required
                for x in range(0,nFiber):
                    
                    if x % 10 == 0:
                        print(str(nFiber-x) + ' Fibers Remaining')
                    
                    # Generate fiber
                    genFibers = self.generateFibers(vertices,tets,coord1,\
                        coord2,coord3,coord4,maxIter,lFiber,maxC,maxAggD,\
                        fiberOrientation,orientationStrength,self.triangles,\
                        cutFiber)
                    self.p1Fibers[x,:] = self.p1Fiber
                    self.p2Fibers[x,:] = self.p2Fiber
                    self.orienFibers[x,:] = self.orienFiber
                    self.fiberLengths[x,:] = self.fiberLength

                fiberTime = round(time.time() - fiberStartTime,2)

                print(str(nFiber) + ' fibers placed in ' + \
                    str(fiberTime) + ' seconds' )
                print('----------------------------------------------------')

        tetTessTimeStart = time.time()

        print('Beginning tetrahedralization.')

        tetGen = self.tetrahedralization(self.vertices2D,\
            self.triangles2D,geoName,geoFile,verbose)

        [allNodes,allTets] = self.readTetgen(Path('meshes/' + geoName + '/' + \
            geoName+'.node'),Path('meshes/' + geoName + '/' + geoName + '.ele'))

        print('Beginning tesselation.')

        [tetFacets,facetCenters,facetAreas,facetNormals,tetn1,tetn2,tetPoints,allDiameters] = \
            self.tesselation(allNodes,allTets,aggDiameterList,minAggD,\
            geoName)

        tetTessTime = round(time.time() - tetTessTimeStart,2)   

        if edgeElements in ['on','On','Y','y','Yes','yes']:

            # If edge elements are turned on, perform edge computations
            edgeData = self.edgeElementData(surfaceExtLength,allNodes,allTets,tetPoints,maxAggD,\
                vertices,tets,coord1,coord2,coord3,coord4,maxC)

        else:

            edgeData = 0

        print('----------------------------------------------------')


        # Leave commented out; section used for study / debugging
        #for materialRule in range(1,11):

        # Leave commented out; section used for study / debugging
        #for materialRule in range(1,11):


        if multiMaterial in ['on','On','Y','y','Yes','yes']:

            edgeMaterialList = 0

        elif cementStructure in ['on','On','Y','y','Yes','yes']:

            edgeMaterialList = self.edgeNodeMaterial(AllHydrationProduct,allNodes,minC,maxC,cementStructureResolution,materialList)
            #edgeMaterialList = 0

        else:

            edgeMaterialList = 0


        if (multiMaterial not in ['on','On','Y','y','Yes','yes']) and (cementStructure not in ['on','On','Y','y','Yes','yes']):

            materialRule = 0


        writeTimeStart = time.time()


        # Adjust values by length scale factor
        if outputUnits in ['nm','NM','nano','nanometer','nanos','nanometers']:

            factor = 1000000.0

        if outputUnits in ['um','UM','micro','micrometer','micros','micrometers']:

            factor = 1000.0

        if outputUnits in ['cm','CM','centi','centimeter','centis','centimeters']:

            factor = 0.1

        if outputUnits in ['dm','DM','deci','decimeter','decis','decimeters']:

            factor = 0.01

        if outputUnits in ['m','M','meter','meters']:

            factor = 0.001

        if outputUnits not in ['mm','MM','milli','millimeter','millis','millimeters','default']:
        
            allNodes = allNodes*factor
            self.nodes = self.nodes*factor
            #cells = [[[x*factor for x in y] for y in z] for z in cells]   
            tetFacets = tetFacets*factor
            facetCenters = facetCenters*factor
            facetAreas = facetAreas*factor**2
            allDiameters = allDiameters*factor
            vertices = vertices*factor
            
            if edgeElements in ['on','On','Y','y','Yes','yes']:
                edgeData[0,0:6] = edgeData[0,0:6]*factor
                edgeData[:,6] = edgeData[:,6]*factor**2
                edgeData[:,10:12] = edgeData[:,10:12]*factor
                edgeData[:,12:14] = edgeData[:,12:14]*factor**3
                edgeData[:,14] = edgeData[:,14]*factor
            
            self.vertices2D = self.vertices2D*factor
            aggDiameterList = aggDiameterList*factor
            minAggD = minAggD*factor
            maxAggD = maxAggD*factor
            
            if fibers in ['on','On','Y','y','Yes','yes']:
                self.p1Fibers = self.p1Fibers*factor
                self.p2Fibers = self.p2Fibers*factor
                dFiber = dFiber*factor
                self.fiberLengths = self.fiberLengths*factor


            # Edit mesh VTK file created by Tetgen
            with open(Path('meshes/' + geoName + '/' + geoName + '-para-mesh.000.vtk')) as myFile:                                      
              for identifier, line in enumerate(myFile, 1):                             
                  if 'POINTS' in line:                                                
                      linePOINTS = int(identifier)                                   
                  if 'CELLS' in line:                                                   
                      lineCELLS = int(identifier)                                        

            points = np.loadtxt(Path('meshes/' + geoName + '/' + geoName + '-para-mesh.000.vtk'), usecols=range(3), \
                skiprows=linePOINTS, max_rows=lineCELLS-linePOINTS-2)                                   

            points = points*factor


            with open(Path('meshes/' + geoName + '/' + geoName + '-para-mesh.000.vtk'), "r") as f:
                contents = f.readlines()

            for x in range(0,len(points)):
                contents.insert(linePOINTS+x, str(points[x,0]) + ' ' \
                    + str(points[x,1]) + ' '  + str(points[x,2]) + '\n')

            del contents[linePOINTS+x+1:lineCELLS-linePOINTS+x+4]

            with open(Path('meshes/' + geoName + '/' + geoName + '-para-mesh.000.vtk'), "w") as f:
                contents = "".join(contents)
                f.write(contents)

        print('Generating Facet Data information.')

        [dataList,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial] = self.facetData(\
            allNodes,allTets,tetFacets,facetCenters,facetAreas,facetNormals,tetn1,\
            tetn2,materialList,materialRule,multiMaterial,cementStructure,edgeMaterialList)

        print('Writing External Facet Data file.')

        # Create file of external triangle facets for plotting of cells
        externalFacetsFile = self.externalFacetFile(dataList,vertices,self.triangles,geoName)

        # Initialize counter for number of facet materials switched
        matSwitched = 0

        # Calculate volume associated with each material (and adjust for high-order material rules)
        if multiMaterial in ['on','On','Y','y','Yes','yes']:

            [itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
                aggVolFracAct] = self.materialVolCheck(geoName,subtetVol,facetMaterial,\
                agg,itz,binder)
            
            if materialRule > 9:

                sortedData1 = self.materialSorting(facetMaterial,facetVol1,facetVol2,\
                    particleMaterial,subtetVol)

                i = 0
                
                while (abs(itzVolFracSim-itzVolFracAct) > 0.02 or abs(itzVolFracSim-itzVolFracAct) == 0.00) and \
                    abs(binderVolFracSim-binderVolFracAct) > 0.02 and \
                    abs(aggVolFracSim-aggVolFracAct) > 0.02 and i < len(sortedData1):

                    # Skip refinement for facets with same-material particles
                    if sortedData1[i,3] != sortedData1[i,4]:
                        
                        # Refine material assignment based on volume fractions
                        sortedData = self.materialRefinement(sortedData1,
                            itzVolFracSim,binderVolFracSim,aggVolFracSim,\
                            itzVolFracAct,binderVolFracAct,aggVolFracAct,i)

                        # Recalculate and update volume fractions
                        [itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
                            aggVolFracAct] = self.materialVolCheck(geoName,sortedData[:,5],sortedData[:,2],\
                            agg,itz,binder)     

                        sortedData1 = sortedData
                        matSwitched = matSwitched+1
                        #print(str(i) + ' / ' + str(len(sortedData)) + ' | ' + str(abs(itzVolFracSim-itzVolFracAct)) + ' | ' + str(abs(binderVolFracSim-binderVolFracAct)) + ' | ' + str(abs(aggVolFracSim-aggVolFracAct)))

                    else:
                        pass

                    i = i+1
                    

                [dataList,facetMaterial] = self.reformDataList(allTets,dataList,sortedData1)

            [itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
                aggVolFracAct] = self.materialVolCheck(geoName,subtetVol,facetMaterial,\
                agg,itz,binder)

            matPlot = self.materialVolPlot(geoName,itzVolFracSim,binderVolFracSim,aggVolFracSim,\
                itzVolFracAct,binderVolFracAct,aggVolFracAct)

            PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
                PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct\
                 = 0,0,0,0,0,0,0,0,0,0

        # Calculate volume associated with each material (and adjust for high-order material rules)
        elif cementStructure in ['on','On','Y','y','Yes','yes']:

            [PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
                PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct]\
                = self.cementMaterialVolCheck(geoName,subtetVol,facetMaterial,PoresVoxelData,\
                ClinkerVoxelData,CHVoxelData,CSH_LDVoxelData,CSH_HDVoxelData,AllHydrationProduct)
            
            if materialRule > 9:

                sortedData1 = self.materialSorting(facetMaterial,facetVol1,facetVol2,\
                    particleMaterial,subtetVol)

                i = 0
                
                while abs(PoresVolFracSim-PoresVolFracAct) > 0.02 and \
                    abs(ClinkerVolFracSim-ClinkerVolFracAct) > 0.02 and \
                    abs(CHVolFracSim-CHVolFracAct) > 0.02 and \
                    abs(CSH_LDVolFracSim-CSH_LDVolFracAct) > 0.02 and \
                    abs(CSH_HDVolFracSim-CSH_HDVolFracAct) > 0.02  and \
                    i < len(sortedData1):

                    # Skip refinement for facets with same-material particles
                    if sortedData1[i,3] != sortedData1[i,4]:
                        
                        # Refine material assignment based on volume fractions
                        sortedData = self.cementMaterialRefinement(sortedData1,
                            PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
                            PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,i)

    
                        # Recalculate and update volume fractions
                        [PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
                            PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct]\
                            = self.cementMaterialVolCheck(geoName,sortedData[:,5],sortedData[:,2],PoresVoxelData,\
                            ClinkerVoxelData,CHVoxelData,CSH_LDVoxelData,CSH_HDVoxelData,AllHydrationProduct)     

                        sortedData1 = sortedData
                        matSwitched = matSwitched+1
                        #print(str(i) + ' / ' + str(len(sortedData)) + ' | ' + str(abs(itzVolFracSim-itzVolFracAct)) + ' | ' + str(abs(binderVolFracSim-binderVolFracAct)) + ' | ' + str(abs(aggVolFracSim-aggVolFracAct)))

                    else:
                        pass

                    i = i+1

                [dataList,facetMaterial] = self.reformDataList(allTets,dataList,sortedData1)

            [PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
                PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct]\
                = self.cementMaterialVolCheck(geoName,subtetVol,facetMaterial,PoresVoxelData,\
                    ClinkerVoxelData,CHVoxelData,CSH_LDVoxelData,CSH_HDVoxelData,AllHydrationProduct)

            #matPlot = self.cementMaterialVolPlot(geoName,PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,\
            #    CSH_LDVolFracSim,CSH_HDVolFracSim,PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,\
            #    CSH_LDVolFracAct,CSH_HDVolFracAct)

            itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct\
             = 0,0,0,0,0,0

        else:

            itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct,\
                PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
                PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,\
                matSwitched,materialRule = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

        if any(elem in output for elem in ['mars', 'MARS', 'Mars']):
            
            print('Writing MARS Data file.')

            # If Mars files requested, generate MARS File
            mkMarsFile = self.marsFile(geoName,allNodes,allTets)

        if any(elem in output for elem in ['abaqus', 'Abaqus', 'ABAQUS']):

            print('Writing Abaqus Data file.')
            
            # If Abaqus files requested, generate Abaqus .inp File
            mkAbaqusFile = self.abaqusFile(geoName,allNodes,allTets,edgeElements,edgeData,self.vertices2D,self.triangles2D)

        if any(elem in output for elem in ['data', 'Data', 'datas', 'Datas']):

            # If data files requested, generate Cell File
            #mkCellFile = self.cellFile(geoName,nParticles,cells,allNodes,dataType)

            print('Writing Mesh Data file.')

            # If data files requested, generate Mesh File
            mkMeshFile = self.meshFile(geoName,allNodes,allTets)

            print('Writing Facet Data file.')

            # If data files requested, generate Facet File
            mkFacetFile = self.facetFile(geoName,dataList,allTets,dataType)

            print('Writing Particle Data file.')

            # If data files requested, generate Particle Data File
            mkParticleData = self.particleData(allNodes,allTets,aggDiameterList,minAggD,geoName)

            if randomField in ['on','On','Y','y','Yes','yes']:
                
                print('Writing Random Field file.')

                # If randomField file requested, generate required Facet File 
                mkFacetFile = self.facetRandomFieldFile(geoName,dataList,allTets,dataType)

            if fibers in ['on','On','Y','y','Yes','yes']:
                
                # If data files requested, generate Fibers File
                mkFibersFile = self.vtkFibers(self.p1Fibers,self.p2Fibers,\
                    dFiber,self.fiberLengths,self.orienFibers,geoName) 

                [FiberdataList,TotalIntersections,MaxInterPerFacet,TotalTet,TotalFiber,IntersectedFiber,projectedFacet]\
                 = self.facetfiberInteractionData(self.p1Fibers,self.p2Fibers,dFiber,self.fiberLengths,self.orienFibers,\
                    geoName,allTets,allNodes,tetFacets,dataList,tetn1,tetn2,facetNormals,facetCenters)
                mkFiberFacetFile = self.facetfiberInteractionFile(geoName,FiberdataList,TotalIntersections,MaxInterPerFacet,TotalTet,TotalFiber,dataType)

                mkNonIntersectedFibersFile = self.vtkNonIntersectedFibers(self.p1Fibers,self.p2Fibers,\
                    dFiber,self.fiberLengths,self.orienFibers,geoName,IntersectedFiber) 

                # If visuals requested, generate Facet VTK File
                mkprojectedFacetVTKFile = self.projectedFacetVTKFile(geoName,tetFacets,dataList,projectedFacet)

            if edgeElements in ['on','On','Y','y','Yes','yes']:

                # Generate edge element data file
                mkEdgeFile = self.edgeFile(geoName,edgeData)

        if any(elem in output for elem in ['visual', 'Visual', 'visuals', 'Visuals']):

            print('Writing visual files.')

            # If visuals requested, generate Particle VTK File
            mkParticleFile = self.vtkParticles(self.nodes,aggDiameterList,\
                materialList,geoName)

            # If visuals requested, generate Facet VTK File
            mkFacetVTKFile = self.facetVTKFile(geoName,tetFacets,dataList,\
                facetMaterial,multiMaterial,cementStructure)

            if edgeElements in ['on','On','Y','y','Yes','yes']:

                # If visuals requested, generate edge element VTK File
                mkEdgeVTKFile = self.edgeVTKFile(geoName,edgeData)

        else:
            
            # If visuals not requested, delete mesh vtk file
            try:
                os.remove(Path('meshes/' + geoName + '/' + geoName + '-para-mesh.000.vtk'))     
            except:
                pass

        # Make sieve curve file
        #mkSieveFile = self.sieveFile(geoName,multiMaterial,aggDiameterList,maxAggD,minAggD,aggGrainsDiameterList,\
        #itzDiameterList,binderDiameterList,cementStructure,PoresDiameterList,ClinkerDiameterList,CHDiameterList,\
        #CSH_LDDiameterList,CSH_HDDiameterList)

        writeTime = round(time.time() - writeTimeStart,2)

        # Generate Log file after run
        mkLogFile = self.logFile(self.gmshTime,nParticles,placementTime,maxAggD,\
            minAggD,fullerCoef,wcRatio,cementC,volFracAir,q,maxIter,\
            geoName,aggOffset,densityWater,densityCement,allTets,dataType,\
            tetTessTime,writeTime,geoFile,dFiber,lFiber,vFiber,fiberFile,\
            multiMaterial,materialFile,maxGrainD,minGrainD,grainFullerCoef,\
            maxBinderD,minBinderD,binderFullerCoef,maxITZD,minITZD,ITZFullerCoef,output,fibers,\
            itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
            aggVolFracAct,sieveCurveDiameter,sieveCurvePassing,matSwitched,materialRule,\
            cementStructure,cementmaterialFile,maxPoresD,minPoresD,PoresFullerCoef,\
            PoresSieveCurveDiameter,PoresSieveCurvePassing,maxClinkerD,minClinkerD,\
            ClinkerFullerCoef,ClinkerSieveCurveDiameter,ClinkerSieveCurvePassing,\
            maxCHD,minCHD,CHFullerCoef,CHSieveCurveDiameter,CHSieveCurvePassing,\
            maxCSH_LDD,minCSH_LDD,CSH_LDFullerCoef,CSH_LDSieveCurveDiameter,CSH_LDSieveCurvePassing,\
            maxCSH_HDD,minCSH_HDD,CSH_HDFullerCoef,CSH_HDSieveCurveDiameter,CSH_HDSieveCurvePassing,\
            PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
            PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,outputUnits)

        # Remove extra files
        os.remove(Path('meshes/' + geoName + '/' + geoName + '.ele'))
        os.remove(Path('meshes/' + geoName + '/' + geoName + '.node'))
        
        # Copy input file as backup
        #shutil.copyfile('LDPMgen.py',Path( 'meshes/' + geoName + '/LDPMgen.py'))

        # Leave commented out; used for study/debugging
        #matLog = self.materialLogFile(materialFile,itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,\
        #       binderVolFracAct,aggVolFracAct,materialRule,minAggD)


####################################################################################################
## Main Class Functions                                                                           ##
####################################################################################################

    def periodicMesher(self,prismX,prismY,prismZ,minAggD,maxAggD,geoName,geoFile):

        start_time_gmsh = time.time()

        # Write the geo file to control periodic meshing
        with open(Path('temp/' + geoName + '.geo'),\
            "w") as f:                                       
            f.write('SetFactory("OpenCASCADE");\n')
            f.write('Box(1) = {0, 0, 0, ' + str(prismX) + ', ' + str(prismY) + ', ' + str(prismZ) + '};\n')
            f.write('MeshSize {:} = ' + str(maxAggD) + ';\n')
            f.write('MeshSize {1} = ' + str(minAggD) + ';\n')
            f.write('Periodic Surface {2} = {1} Translate {' + str(prismX) + ', 0, 0};\n')
            f.write('Periodic Surface {6} = {5} Translate {0, 0, ' + str(prismZ) + '};\n') 
            f.write('Periodic Surface {4} = {3} Translate {0, ' + str(prismY) + ', 0};\n')  
            f.write('\n')


        # Call the geo file for meshing
        print('----------------------------------------------------')
        print('Starting Gmsh with element size ' + str(minAggD) + ' mm to ' + \
            str(maxAggD) + ' mm.')
        gmshCommand = str(Path('lib/gmsh/gmsh')) + ' ' + str(Path('temp/' \
            + geoName + '.geo')) + ' -3 .mesh -format mesh -v 0'                                    
        os.system(gmshCommand)
        
        # Mesh Input - 2D Mesh Conversion
        gmshCommand2D = str(Path('lib/gmsh/gmsh')) + ' ' \
        + str(Path('temp/' + geoFile)) + ' -2 \
         -algo del2d -o ' + str(Path('temp/' + geoName)) + \
         '2D.stl -format stl -v 0'                                    
        os.system(gmshCommand2D)

        gmshCommand2D = str(Path('lib/gmsh/gmsh')) + ' ' + str(Path('temp/'\
            + geoName)) + '2D.stl -2 \
         -algo del2d -o ' + str(Path('temp/' + geoName)) + '2D.mesh -format mesh -v 0'                                    
        os.system(gmshCommand2D)        

        ParticleGen.gmshTime = round(time.time() - start_time_gmsh,2)
        print('Mesh creation completed in ' + str(ParticleGen.gmshTime) \
            + ' seconds.')
        print('----------------------------------------------------')


    def meshConcreteGeometry(self, geoName,minAggD,maxAggD,geoFile,geoExtension):

        start_time_gmsh = time.time()

        if not os.path.exists(Path('temp')):
            os.mkdir(Path('temp'))  

        if geoExtension in ['step', 'STEP', 'stp', 'STP', 'igs', 'IGS', 'iges',\
            'IGES', 'brep', 'BREP', 'geo', 'GEO']:
            
            # Geometry Input - 3D Meshing
            print('----------------------------------------------------')
            print('Starting Gmsh with element size ' + str(minAggD) + ' mm to ' + \
                str(2*minAggD) + ' mm.')
            gmshCommand = str(Path('lib/gmsh/gmsh')) + ' ' + str(Path('geometry/' \
                + geoFile)) + ' -3 -clmin ' + str(minAggD) + ' -clmax ' + str(2*minAggD) + \
              ' -algo del3d -o ' + str(Path('temp/' + geoName)) + '.mesh -format mesh -v 0'                                    
            os.system(gmshCommand)

            # Geometry Input - 2D Meshing
            gmshCommand2D = str(Path('lib/gmsh/gmsh')) + ' ' + str(Path('geometry/' + geoFile)) + ' -2 \
            -clmin ' + str(minAggD) + ' -clmax ' + str(2*minAggD) + \
              ' -algo del2d -o ' + str(Path('temp/' + geoName)) + '2D.mesh -format mesh -v 0'                                    
            os.system(gmshCommand2D)
            ParticleGen.gmshTime = round(time.time() - start_time_gmsh,2)
            print('Meshing completed in ' + str(ParticleGen.gmshTime) + ' seconds.')
            print('----------------------------------------------------')

        elif geoExtension in ['msh', 'MSH', 'diff', 'DIFF', 'unv', 'UNV', \
            'med', 'MED', 'mmed', 'MMED','mesh', 'MESH', 'p3d', 'P3D', 'vtk', \
            'VTK', 'bdf', 'BDF', 'nas', 'NAS']:

            # Mesh Input - 3D Mesh Conversion
            print('----------------------------------------------------')
            print('Starting Gmsh for mesh conversion.')
            gmshCommand = str(Path('lib/gmsh/gmsh')) + ' ' + str(Path('geometry/' \
                + geoFile)) + ' -3 -algo del3d -o ' + str(Path('temp/' \
                + geoName)) + '.mesh -format mesh -v 0'                                    
            os.system(gmshCommand)
            
            # Mesh Input - 2D Mesh Conversion
            gmshCommand2D = str(Path('lib/gmsh/gmsh')) + ' ' \
            + str(Path('geometry/' + geoFile)) + ' -2 \
             -algo del2d -o ' + str(Path('temp/' + geoName)) + \
             '2D.stl -format stl -v 0'                                    
            os.system(gmshCommand2D)

            gmshCommand2D = str(Path('lib/gmsh/gmsh')) + ' ' + str(Path('temp/'\
                + geoName)) + '2D.stl -2 \
             -algo del2d -o ' + str(Path('temp/' + geoName)) + '2D.mesh -format mesh -v 0'                                    
            os.system(gmshCommand2D)        

            ParticleGen.gmshTime = round(time.time() - start_time_gmsh,2)
            print('Mesh conversion completed in ' + str(ParticleGen.gmshTime) \
                + ' seconds.')
            print('----------------------------------------------------')

        else:
            print('The geometry file is not a compatible CAD or mesh format.')
            print('Please read the README and submit with a new file.')
            print('Now exitting...')
            exit()

    def meshRebarGeometry(self,minAggD,maxAggD,rebar,rebarFile1,dRebar1,\
        rebarFile2,dRebar2,rebarFile3,dRebar3,geoName):

        start_time_gmsh = time.time()

        try:
            rebarName1 = rebarFile1.split(".")[0]
            rebarExtension1 = rebarFile1.split(".")[1]
        except:
            print('Rebar is turned on but no files were provided')
            print('Exitting now...')
            exit()

        try:
            rebarName2 = rebarFile2.split(".")[0]
            rebarExtension2 = rebarFile2.split(".")[1]
        except:
            rebarName2 = ''
            rebarExtension2 = ''

        try:
            rebarName3 = rebarFile3.split(".")[0]
            rebarExtension3 = rebarFile3.split(".")[1]      
        except:
            rebarName3 = ''
            rebarExtension3 = ''

        dRebar = [dRebar1,dRebar2,dRebar3]
        rebarFile = [rebarFile1,rebarFile2,rebarFile3]
        rebarName = [rebarName1,rebarName2,rebarName3]
        rebarExtension = [rebarExtension1,rebarExtension2,rebarExtension3]

        if not os.path.exists(Path('temp')):
            os.mkdir(Path('temp'))  

        for x in range(3):

            if rebarExtension[x] == '':
                pass

            elif rebarExtension[x] in ['step', 'STEP', 'stp', 'STP', 'igs', 'IGS', 'iges',\
                'IGES', 'brep', 'BREP', 'geo', 'GEO']:
                
                # Geometry Input - 2D Meshing
                print('----------------------------------------------------')
                print('Starting Gmsh for rebar with size ' + str(1.5*maxAggD) + ' mm to ' + \
                    str(3*maxAggD) + ' mm.')

                gmshCommand2D = str(Path('lib/gmsh/gmsh')) + ' ' + str(Path('rebar/' + rebarFile[x])) + ' -1 \
                -clmin ' + str(1.5*maxAggD) + ' -clmax ' + str(3*maxAggD) + \
                  ' -o ' + str(Path('meshes/' + geoName + '/' + geoName + '-' + rebarName[x])) + '.mesh -format mesh -v 0'                                    
                
                os.system(gmshCommand2D)
                gmshTime = round(time.time() - start_time_gmsh,2)
                print('Meshing completed in ' + str(gmshTime) + ' seconds.')
                print('----------------------------------------------------')

                rebarMesh = Path('meshes/' + geoName + '/' + geoName + '-' + rebarName[x] + '.mesh')

                with open(rebarMesh) as myFile:                                      
                  for identifier, line in enumerate(myFile, 1):                             
                      if 'Vertices' in line:                                                
                          lineVertices = int(identifier)                                   
                      if 'Edges' in line:                                                   
                          lineEdges = int(identifier)                           
                      if 'End' in line:                                                      
                          lineEnd = int(identifier)                                          

                vertices = np.loadtxt(rebarMesh, usecols=range(3), \
                    skiprows=lineVertices+1, max_rows=lineEdges-lineVertices-2)                                   

                edges = np.loadtxt(rebarMesh, usecols=range(3),\
                    skiprows=lineEdges+1, max_rows=lineEnd-lineEdges-2)   

                mkRebarFile = self.vtkRebar(vertices[edges[:,0].astype(int)-1],vertices[edges[:,1].astype(int)-1],dRebar[x],geoName,x)

            elif rebarExtension[x] in ['msh', 'MSH', 'diff', 'DIFF', 'unv', 'UNV', \
                'med', 'MED', 'mmed', 'MMED','mesh', 'MESH', 'p3d', 'P3D', 'vtk', \
                'VTK', 'bdf', 'BDF', 'nas', 'NAS']:

                # Mesh Input - 2D Mesh Conversion
                print('----------------------------------------------------')
                print('Starting Gmsh for rebar mesh conversion.')

                gmshCommand2D = str(Path('lib/gmsh/gmsh')) + ' ' \
                + str(Path('rebar/' + rebarFile[x])) + ' -1 \
                 -o ' + str(Path('temp/' + rebarName[x])) + \
                 '.stl -format stl -v 0'                                    
                os.system(gmshCommand2D)
                
                gmshCommand2D = str(Path('lib/gmsh/gmsh')) + ' ' + str(Path('temp/'\
                    + rebarName[x])) + '.stl -1 \
                 -o ' + str(Path('meshes/' + geoName + '/' + geoName + '-' + rebarName[x])) + '.mesh -format mesh -v 0'                                    
                os.system(gmshCommand2D)    
                os.remove(Path('temp/' + rebarName[x])) 

                gmshTime = round(time.time() - start_time_gmsh,2)
                print('Mesh conversion completed in ' + str(gmshTime) \
                    + ' seconds.')
                print('----------------------------------------------------')

            else:
                print('The geometry file is not a compatible CAD or mesh format.')
                print('Please read the README and submit with a new file.')
                print('Now exitting...')
                exit()

    def readConcreteMesh(self, concreteMeshFile, concreteMeshFile2D):

        with open(concreteMeshFile) as myFile:                                      
          for identifier, line in enumerate(myFile, 1):                             
              if 'Vertices' in line:                                                
                  lineVertices = int(identifier)                                   
              if 'Edges' in line:                                                   
                  lineEdges = int(identifier)   
              if 'Triangles' in line:                                                   
                  lineTriangles = int(identifier)                                      
              if 'Tetrahedra' in line:                                             
                  lineTetrahedra = int(identifier)                            
              if 'End' in line:                                                      
                  lineEnd = int(identifier)                                          

        vertices = np.loadtxt(concreteMeshFile, usecols=range(3), \
            skiprows=lineVertices+1, max_rows=lineEdges-lineVertices-2)                                   

        ParticleGen.triangles = np.loadtxt(concreteMeshFile, usecols=range(3),\
            skiprows=lineTriangles+1, max_rows=lineTetrahedra-lineTriangles-2)   

        ParticleGen.tets = np.loadtxt(concreteMeshFile, usecols=range(4), \
            skiprows=lineTetrahedra+1, max_rows=lineEnd-lineTetrahedra-2)       

        lineEdges = 0

        with open(concreteMeshFile2D) as myFile:                                      
          for identifier, line in enumerate(myFile, 1):                             
              if 'Vertices' in line:                                                
                  lineVertices = int(identifier)                                   
              if 'Edges' in line:                                                   
                  lineEdges = int(identifier) 
              if 'Triangles' in line:                                                   
                  lineTriangles = int(identifier)                                                              
              if 'End' in line:                                                      
                  lineEnd = int(identifier)                                          

        if lineEdges == 0:
            ParticleGen.vertices2D = np.loadtxt(concreteMeshFile2D, \
                usecols=range(4), skiprows=lineVertices+1, \
                max_rows=lineTriangles-lineVertices-2)      

        else: 
            ParticleGen.vertices2D = np.loadtxt(concreteMeshFile2D, \
                usecols=range(4), skiprows=lineVertices+1, \
                max_rows=lineEdges-lineVertices-2)      

            ParticleGen.edges2D = np.loadtxt(concreteMeshFile2D, \
                usecols=range(3), skiprows=lineEdges+1, \
                max_rows=lineTriangles-lineEdges-2)                  

        ParticleGen.triangles2D = np.loadtxt(concreteMeshFile2D, \
            usecols=range(4), skiprows=lineTriangles+1, \
            max_rows=lineEnd-lineTriangles-2)   

        return vertices

    # Read the voxelated microstructure file and store
    def readMultiMaterial(self, materialFile):

        microstructureFile = Path('microstructure/' + materialFile)

        with open(microstructureFile) as f:
            version = str(f.readlines(1))
            X_Size = str(f.readlines(2))
            Y_Size = str(f.readlines(3))
            Z_Size = str(f.readlines(4))
            resolution = str(f.readlines(5))

        X_Size = X_Size[10:]
        Y_Size = Y_Size[10:]
        Z_Size = Z_Size[10:]
        resolution = resolution[20:]

        X_Size = float(X_Size[:-4])
        Y_Size = float(Y_Size[:-4])
        Z_Size = float(Z_Size[:-4])
        resolution = float(resolution[:-4])

        voxels = np.loadtxt(microstructureFile, skiprows=5)

        return X_Size, Y_Size, Z_Size, resolution, voxels

    # Check if the voxelated microstructure is larger than the geometry
    def checkMicrostructureSize(self,X_Size,Y_Size,Z_Size,resolution,minC,maxC):
        microstructureSize = np.array((X_Size,Y_Size,Z_Size))*resolution
        geoSize = maxC-minC

        if ((microstructureSize-geoSize)>=0).all():
            return True
        else:
            print('The microstructure is smaller than the geometry.')
            print('Now exitting...')
            exit()

    def sortMicrostructure(self,voxels):

        # Number all voxels
        voxelNumbering = np.arange(len(voxels))+1

        # Keep only aggregate voxels
        aggFull = (voxels>3)*voxelNumbering
        agg = aggFull[aggFull != 0]
        
        # Keep only ITZ voxels
        itzFull = (voxels==2)*voxelNumbering
        itz = itzFull[itzFull != 0]

        # Keep only binder voxels
        binderFull = (voxels==3)*voxelNumbering
        binder = binderFull[binderFull != 0]

        return agg,itz,binder

    # Read the voxelated cementStructure file and store
    def readCementStructure(self, cementmaterialFile, minC, cementStructureResolution):

        cementStructureFile = Path('cementstructure/' + cementmaterialFile)
        #if VoxelNumbering>100 then VoxelNumbering must be a multiple of 100
        CementStructureVoxels = np.loadtxt(cementStructureFile, usecols=range(1))
        VoxelNumbering = np.arange(len(CementStructureVoxels)).reshape((-1,1))
        CementStructure_X_Size = round((len(VoxelNumbering))**(1/3))
        CementStructure_Y_Size = round((len(VoxelNumbering))**(1/3))
        CementStructure_Z_Size = round((len(VoxelNumbering))**(1/3))

        VoxelCenterX = np.zeros(len(VoxelNumbering))
        VoxelCenterY = np.zeros(len(VoxelNumbering))
        VoxelCenterZ = np.zeros(len(VoxelNumbering))
        XYZID = []
        for k in range(len(VoxelNumbering)):
            # Find voxel coordinates of voxel k 
            xVoxel = np.floor(k/(CementStructure_Z_Size*CementStructure_Y_Size))+1
            yVoxel = np.floor((k-(xVoxel-1)*(CementStructure_Z_Size*CementStructure_Y_Size))/CementStructure_Z_Size)+1
            zVoxel = k-(xVoxel-1)*(CementStructure_Z_Size*CementStructure_Y_Size)-((yVoxel-1)*CementStructure_Z_Size)+1

            VoxelCenterX[k]=int(xVoxel-1)*cementStructureResolution+cementStructureResolution/2+minC[0]
            VoxelCenterY[k]=int(yVoxel-1)*cementStructureResolution+cementStructureResolution/2+minC[1]
            VoxelCenterZ[k]=int(zVoxel-1)*cementStructureResolution+cementStructureResolution/2+minC[2]

            DataRow = np.array([int(k+1), VoxelCenterX[k], VoxelCenterY[k], VoxelCenterZ[k], CementStructureVoxels[k]])
            XYZID.append(DataRow)
        XYZID=np.array(XYZID).reshape(-1,5)
        print(CementStructure_X_Size, CementStructure_Y_Size, CementStructure_Z_Size)

        return CementStructure_X_Size, CementStructure_Y_Size, CementStructure_Z_Size, XYZID, CementStructureVoxels

    # Check if the voxelated cementStructure is larger than the geometry
    def checkCementStructureSize(self,CementStructure_X_Size,CementStructure_Y_Size,CementStructure_Z_Size,cementStructureResolution,minC,maxC):
        CementStructureSize = np.array((CementStructure_X_Size,CementStructure_Y_Size,CementStructure_Z_Size))*cementStructureResolution
        geoSize = maxC-minC

        if ((CementStructureSize-geoSize)>=0).all():
            return True
        else:
            print('The cementStructure is smaller than the geometry.')
            print('Now exitting...')
            exit()

    def sortCementStructure(self,CementStructureVoxels,XYZID,wcRatio,Alphac,SaturatedCSHDensity,Mean_CSH_HD,Mean_CSH_LD,SDev_CSH_HD,SDev_CSH_LD,densityWater):

        # Number all voxels
        VoxelNumbering = np.arange(len(CementStructureVoxels))+1
        PoresVoxelData = []
        ClinkerVoxelData = []
        OtherVoxelData = []
        CHVoxelData = []
        CSHVoxelData = []
        CSH_LDVoxelData = []
        CSH_HDVoxelData = []        

        for i in range(len(VoxelNumbering)):

        # Coordinates of Capillary Pores voxels with properties to Capillary Pores voxels 
            if XYZID[i,4] == 0 or XYZID[i,4] == 45:
                CAPPORData = np.array([XYZID[i,0], XYZID[i,1], XYZID[i,2], XYZID[i,3], int(1)])
                PoresVoxelData.append(CAPPORData)

        # Coordinates of Clinker voxels with properties to Clinker voxels 
            elif XYZID[i,4] >= 1 and XYZID[i,4] <= 12:
                ClinkerData = np.array([XYZID[i,0], XYZID[i,1], XYZID[i,2], XYZID[i,3], int(2)])
                ClinkerVoxelData.append(ClinkerData)

        # Coordinates of CSH voxels
            elif XYZID[i,4] == 14 or  XYZID[i,4] == 20 or XYZID[i,4] == 21 or XYZID[i,4] == 30:
                CSHData = np.array([XYZID[i,0], XYZID[i,1], XYZID[i,2], XYZID[i,3], int(4)])
                CSHVoxelData.append(CSHData)

        # Coordinates of CH voxels with properties to CH voxels 
            elif XYZID[i,4] == 13 or  XYZID[i,4] == 31:
                CHData = np.array([XYZID[i,0], XYZID[i,1], XYZID[i,2], XYZID[i,3], int(3)])
                CHVoxelData.append(CHData)
        # Coordinates of CH voxels with properties to CH voxels 
            else:                 
                OtherData = np.array([XYZID[i,0], XYZID[i,1], XYZID[i,2], XYZID[i,3], int(4)])
                OtherVoxelData.append(OtherData)

        PoresVoxelData = np.array(PoresVoxelData).reshape(-1,5)
        ClinkerVoxelData = np.array(ClinkerVoxelData).reshape(-1,5)
        CHVoxelData = np.array(CHVoxelData).reshape(-1,5)
        CSHVoxelData = np.array(CSHVoxelData).reshape(-1,5)
        OtherVoxelData = np.array(OtherVoxelData).reshape(-1,5)
        NumberofCSHVoxel = len(CSHVoxelData)
    
        M_r = 3.017 * wcRatio * Alphac - 1.347 * Alphac + 0.538 
        if M_r < 0.0:
            M_r = 0.01
        #print(M_r)

        Density_LD = SaturatedCSHDensity*(1-(1-Mean_CSH_LD)*(1-densityWater/SaturatedCSHDensity))
        Density_HD = SaturatedCSHDensity*(1-(1-Mean_CSH_HD)*(1-densityWater/SaturatedCSHDensity))

        VolFracLD_CSH = (M_r*(SaturatedCSHDensity*Mean_CSH_HD+densityWater*(1-Mean_CSH_HD)))/\
                            (M_r*(Mean_CSH_HD-Mean_CSH_LD)+(SaturatedCSHDensity-densityWater)+\
                                SaturatedCSHDensity*Mean_CSH_LD+densityWater*(1-Mean_CSH_LD))

        VolFracHD_CSH = 1.0 - VolFracLD_CSH
        #print(VolFracHD_CSH)

        # Generate random numbers
        randomNCSH = np.random.rand(NumberofCSHVoxel)
        VolumnFraction = np.array([VolFracLD_CSH, VolFracHD_CSH])
        MeanValue = np.array([Mean_CSH_LD, Mean_CSH_HD])
        StandardDeviation = np.array([SDev_CSH_LD, SDev_CSH_HD])
        eta_min = 0.5
        eta_max = 1.0
        PakingDensity = np.zeros(len(randomNCSH))
        Density = np.zeros(len(randomNCSH))
        Pores = np.zeros(len(randomNCSH))
        CSHID = np.zeros(len(randomNCSH))
        j = 0
        S1 = StandardDeviation[j]
        M1 = MeanValue[j]
        f1 = VolumnFraction[j]
        S2 = StandardDeviation[j+1]
        M2 = MeanValue[j+1]
        f2 = VolumnFraction[j+1]

        def integrandLD(a):
            return np.exp(-1*np.power((a-M1),2)/(2*np.power(S1,2)))
        def integrandHD(a):
            return np.exp(-1*np.power((a-M2),2)/(2*np.power(S2,2)))

        for k in range(len(randomNCSH)):
            def All(eta):
                return abs(quad(integrandLD, np.NINF, eta)[0]*f1*(1/(S1*math.sqrt(2*np.pi)))+\
                quad(integrandHD, np.NINF, eta)[0]*f2*(1/(S2*math.sqrt(2*np.pi)))-randomNCSH[k])

            PakingDensity[k] = scipy.optimize.fminbound(All, eta_min, eta_max)
            Pores[k] = 1-PakingDensity[k]
            Density[k] = SaturatedCSHDensity*(1-Pores[k]*(1-densityWater/SaturatedCSHDensity))
            #print(PakingDensity[k])
            #print((M2-S2))

        # Coordinates of CSH_LD voxels with properties to Low Density CSH voxels 
            if PakingDensity[k] <= (M2-S2):
                CSHID[k] = int(4)
                LD_CSHData = np.array([CSHVoxelData[k,0], CSHVoxelData[k,1], CSHVoxelData[k,2], CSHVoxelData[k,3], CSHID[k]])
                CSH_LDVoxelData.append(LD_CSHData)
        # Coordinates of CSH_HD voxels with properties to Height Density CSH voxels 
            else:
                CSHID[k] = int(5)
                #print(CSHID[k])
                HD_CSHData = np.array([CSHVoxelData[k,0], CSHVoxelData[k,1], CSHVoxelData[k,2], CSHVoxelData[k,3], CSHID[k]])
                CSH_HDVoxelData.append(HD_CSHData)

        CSH_LDVoxelData = np.array(CSH_LDVoxelData).reshape(-1,5)
        CSH_LDVoxelData = np.vstack((CSH_LDVoxelData , OtherVoxelData))
        CSH_HDVoxelData = np.array(CSH_HDVoxelData).reshape(-1,5)

        AllHydrationProduct = np.concatenate((PoresVoxelData, ClinkerVoxelData, CHVoxelData, CSH_LDVoxelData, CSH_HDVoxelData), axis=0)

        return PoresVoxelData, ClinkerVoxelData, CHVoxelData, CSH_LDVoxelData, CSH_HDVoxelData, AllHydrationProduct

    # Give material to the edge nodes
    def edgeNodeMaterial(self,AllHydrationProduct,allNodes,minC,maxC,cementStructureResolution,materialList):
        edgeVoxels = []
        edgeNodes = []
        edgeMaterialList =[]
        eps = 10**(-8)

        List = int(abs(len(allNodes)-len(materialList)))
        Edgelist = int(len(edgeNodes))
        #print(len(materialList))

        for i in range(0,len(AllHydrationProduct)):

            minedgeX = AllHydrationProduct[i,1]-cementStructureResolution/2
            maxedgeX = AllHydrationProduct[i,1]+cementStructureResolution/2
            minedgeY = AllHydrationProduct[i,2]-cementStructureResolution/2
            maxedgeY = AllHydrationProduct[i,2]+cementStructureResolution/2
            minedgeZ = AllHydrationProduct[i,3]-cementStructureResolution/2
            maxedgeZ = AllHydrationProduct[i,3]+cementStructureResolution/2
            MaterialID = AllHydrationProduct[i,4]

            edgeV = np.array([minedgeX, maxedgeX, minedgeY, maxedgeY, minedgeZ, maxedgeZ, MaterialID])
            edgeVoxels.append(edgeV)

        edgeVoxels=np.array(edgeVoxels).reshape(-1,7)
        #print(edgeVoxels.shape)

        i=0
        while (Edgelist < List):
            tempEdgeNodes = []
            newAllNodes = []

            #Convexhull outer
            ConvexhullVertices = np.array(ConvexHull(allNodes).vertices).reshape(-1,1)
            print(allNodes.shape)

            for j in range(0,len(allNodes)):

                Xnode = allNodes[j,0]
                Ynode = allNodes[j,1]
                Znode = allNodes[j,2]
                edgeN = np.array([Xnode, Ynode, Znode])

                if j in ConvexhullVertices :
                    tempEdgeNodes.append(edgeN)
                else:
                    newAllNodes.append(edgeN)

            tempEdgeNodes=np.array(tempEdgeNodes).reshape(-1,3)
            #print(tempEdgeNodes.shape)
            newAllNodes=np.array(newAllNodes).reshape(-1,3)
            #print(newAllNodes.shape)

            if i==0:
                edgeNodes=tempEdgeNodes

            else:
                edgeNodes = np.concatenate((tempEdgeNodes,edgeNodes))

            Edgelist = int(len(edgeNodes))

            allNodes = newAllNodes
            i=i+1

        #print(edgeNodes.shape)

        # Remove extra 
        # while (Edgelist > List):
        #     edgeNodes = np.delete(edgeNodes, -1, axis=0)

        # if (Edgelist > List):
        #     edgeNodes = np.delete(edgeNodes, slice(int(List+1), -1), axis=0)

        if (Edgelist > List):
            edgeNodes = edgeNodes[:List]
        #print(edgeNodes.shape)

        for K in range(0,len(edgeNodes)):

                for L in range(0,len(edgeVoxels)):

                    if  ((edgeVoxels[L,0]<=(edgeNodes[K,0]+eps) and (edgeNodes[K,0]-eps)<=edgeVoxels[L,1]) and\
                        (edgeVoxels[L,2]<=(edgeNodes[K,1]+eps) and (edgeNodes[K,1]-eps)<=edgeVoxels[L,3]) and\
                        (edgeVoxels[L,4]<=(edgeNodes[K,2]+eps) and (edgeNodes[K,2]-eps)<=edgeVoxels[L,5])):

                        edgeMaterial = np.array([edgeVoxels[L,6]])
                        edgeMaterialList.append(edgeMaterial)
                        break
                    else: continue


        edgeMaterialList = np.array(edgeMaterialList).reshape(-1)
        #print(edgeMaterialList)
        print(edgeMaterialList.shape)

        return edgeMaterialList
    # Pulls the coordinates of each external triangle in the mesh 
    def formTriangles(self, vertices, triangles):

        triangles = triangles.astype(int)
        
        coord1 = vertices[triangles[:,0]-1]
        coord2 = vertices[triangles[:,1]-1]
        coord3 = vertices[triangles[:,2]-1] 

        trianglePoints = np.concatenate((coord1,coord2,coord3))

        return trianglePoints


    def meshExtents(self,vertices):

        minC = np.amin(vertices, axis=0)
        maxC = np.amax(vertices, axis=0)

        return minC, maxC


    def meshVolume(self,vertices,tets):

        tets = tets.astype(int)
        
        coord1 = vertices[tets[:,0]-1]
        coord2 = vertices[tets[:,1]-1]
        coord3 = vertices[tets[:,2]-1]
        coord4 = vertices[tets[:,3]-1]

        tetVolume = abs(np.vdot(np.transpose(coord1-coord4),\
            np.transpose(np.cross((coord2-coord4),(coord3-coord4)))))/6

        return tetVolume


    # Shift total sieve curve to coarse sieve curve
    def sieveCurve(self,minAggD,maxAggD,sieveCurveDiameter,sieveCurvePassing):
        
        nblines = len(sieveCurveDiameter)

        for i in range(0,nblines):
            if minAggD >= sieveCurveDiameter[i] and minAggD < sieveCurveDiameter[i+1]:      
                minRange = i

        for i in range(0,nblines):
            if maxAggD > sieveCurveDiameter[i] and maxAggD <= sieveCurveDiameter[i+1]:              
                maxRange = i

        w_min = sieveCurvePassing[minRange]+(sieveCurvePassing[minRange+1]-sieveCurvePassing[minRange])/\
            (sieveCurveDiameter[minRange+1]-sieveCurveDiameter[minRange])*(minAggD-sieveCurveDiameter[minRange])
        w_max = sieveCurvePassing[maxRange]+(sieveCurvePassing[maxRange+1]-sieveCurvePassing[maxRange])/\
            (sieveCurveDiameter[maxRange+1]-sieveCurveDiameter[maxRange])*(maxAggD-sieveCurveDiameter[maxRange])

        NewSet = maxRange-minRange+1

        newSieveCurveD = [0]*100
        newSieveCurveP = [0]*100

        for i in range(0,NewSet+1):
            if i == 0:
                newSieveCurveD[i] = minAggD
                newSieveCurveP[i] = 0
            elif NewSet > 1 and i > 0 and i < NewSet:
                newSieveCurveD[i] = sieveCurveDiameter[minRange + i]
                newSieveCurveP[i] = (sieveCurvePassing[minRange + i] - w_min)/(w_max - w_min)
            elif NewSet == i:
                newSieveCurveD[i] = maxAggD
                newSieveCurveP[i] = 1

        newSieveCurveD = np.trim_zeros(newSieveCurveD,trim='b')
        newSieveCurveP = np.trim_zeros(newSieveCurveP,trim='b')

        return newSieveCurveD,newSieveCurveP,NewSet,w_min,w_max


    def aggVolume(self,tetVolume,wcRatio,cementC,volFracAir,fullerCoef,\
        densityCement,densityWater,minAggD,maxAggD,newSieveCurveD,newSieveCurveP,NewSet,w_min,w_max):



        volFracCement = cementC/densityCement
        volFracWater = wcRatio*cementC/densityWater
        volFracAgg = (1-volFracCement-volFracWater-volFracAir)


        if newSieveCurveD == 0:

            volFracAggSim = (1-(minAggD/maxAggD)**fullerCoef)*volFracAgg

            cdf = 0
            cdf1 = 0
            kappa_i = 0

        else:

            # Calculate Kappa
            A = 0 

            for i in range(0,NewSet):
                B = (newSieveCurveP[i+1] - newSieveCurveP[i])/(newSieveCurveD[i+1] - newSieveCurveD[i])\
                    *(1/newSieveCurveD[i]/newSieveCurveD[i] - 1/newSieveCurveD[i+1]/newSieveCurveD[i+1])
                A = A+B

            kappa = 2/A
            kappa_i = [0]*100

            for i in range(0,NewSet):
                kappa_i[i] = (kappa * (newSieveCurveP[i+1] - newSieveCurveP[i]) / (newSieveCurveD[i+1] - newSieveCurveD[i]))
            kappa_i = np.trim_zeros(kappa_i,trim='b')

            # Calculate Volume Fraction
            w_sim = w_max-w_min
            volFracAggSim = w_sim*volFracAgg

            # Get CDF
            cdf = [0]*100
            cdf1 = [0]*100
            for i in range(0,NewSet):
                cdf1[i] = kappa_i[i] * (1./(newSieveCurveD[i]**2) - 1/(newSieveCurveD[i+1]**2))/2
                if i > 0:
                    cdf[i] = cdf1[i-1] + cdf[i-1]
            cdf[NewSet] = cdf1[NewSet-1] + cdf[NewSet-1]

            cdf = np.trim_zeros(cdf,trim='b')
            cdf1 = np.trim_zeros(cdf1,trim='b')

        aggVolTotal = volFracAggSim*tetVolume

        return aggVolTotal,cdf,cdf1,kappa_i


    def aggList(self,aggVolTotal,minAggD,maxAggD,q,newSieveCurveD,cdf,kappa_i,NewSet):

        smallAggVolume = 4/3*math.pi*(minAggD/2)**3
        maxAggNum = np.ceil(aggVolTotal/smallAggVolume)
        aggDiameter = np.zeros(int(maxAggNum))
        aggVol = np.zeros(int(maxAggNum))

        i = 0

        # Fuller Coefficient Case
        if newSieveCurveD == 0:

            # Count particles individually
            if len(aggDiameter) <= 100:
                while sum(aggVol) < aggVolTotal:
                    aggDiameter[i] = minAggD*(1-np.random.rand(1)*(1-minAggD**q\
                        /maxAggD**q))**(-1/q)
                    aggVol[i] = 4/3*math.pi*(aggDiameter[i]/2)**3
                    i = i+1         

            # Count particles in bunches of 100
            elif len(aggDiameter) <= 1000:
                while sum(aggVol) < aggVolTotal:
                    aggDiameter[i:i+100] = minAggD*(1-np.random.rand(100)*\
                        (1-minAggD**q/maxAggD**q))**(-1/q)
                    aggVol[i:i+100] = 4/3*math.pi*(aggDiameter[i:i+100]/2)**3
                    i = i+100

            # Count particles in bunches of 1000
            else:
                while sum(aggVol) < aggVolTotal:
                    aggDiameter[i:i+1000] = minAggD*(1-np.random.rand(1000)*\
                        (1-minAggD**q/maxAggD**q))**(-1/q)
                    aggVol[i:i+1000] = 4/3*math.pi*(aggDiameter[i:i+1000]/2)**3
                    i = i+1000

        # Sieve Curve Case
        else:

            while sum(aggVol) < aggVolTotal:

                F = np.random.rand(1)

                for x in range(0,NewSet):

                    if (F >= cdf[x] and F < cdf[x+1]) :
                        aggDiameter[i] = ((newSieveCurveD[x]**2*kappa_i[x])/(kappa_i[x]-2*(F-cdf[x])*newSieveCurveD[x]**2))**0.5
                        aggVol[i] = 4/3*math.pi*(aggDiameter[i]/2)**3
                        i = i+1         

        # Remove trailing zeros from arrays
        aggDiameter = np.trim_zeros(aggDiameter)
        aggVol = np.trim_zeros(aggVol)

        # Remove accidental extra placed particles
        while sum(aggVol) > aggVolTotal:
            aggDiameter = np.delete(aggDiameter, -1)
            aggVol = np.delete(aggVol, -1)

        # Sort particle diameters large-to-small
        aggDiameterList = np.sort(aggDiameter)[::-1]

        return maxAggNum,aggDiameterList


    def generateParticle(self,x,trianglePoints,aggDiameter,maxAggNum,minC,maxC,\
        vertices,tets,coord1,coord2,coord3,coord4,newMaxIter,maxIter,minAggD,maxAggD,\
        aggOffset,verbose,aggDiameterList,maxEdgeLength):
        
        # Generate random numbers to use in generation
        randomN = np.random.rand(newMaxIter*3)
        i=0
        ntet = len(tets)
        # Generate random nodal location
        while True:
            i=i+3

            if i/3 >= newMaxIter:
                i = 0
                newMaxIter = newMaxIter*2
                randomN = np.random.rand(newMaxIter*3)

            if newMaxIter >= maxIter:
                print("This particle has exceeeded the %r specified maximum iterations allowed." % (maxIter))
                print('Now exitting...')
                exit()

            # Random point selection in random tet prism container
            tetVerts = np.vstack((vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),0]-1),:],\
                vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),1]-1),:],\
                vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),2]-1),:],\
                vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),3]-1),:]))

            tetMin = np.amin(tetVerts, axis=0)
            tetMax = np.amax(tetVerts, axis=0)

            node = np.array([randomN[i]*(tetMax[0]-tetMin[0])+tetMin[0],\
                randomN[i+1]*(tetMax[1]-tetMin[1])+tetMin[1],randomN[i+2]\
                *(tetMax[2]-tetMin[2])+tetMin[2]]).T
            node = node[np.newaxis,:]           


            # Obtain extents for floating bin
            binMin = np.array(([node[0,0]-aggDiameter/2-maxAggD/2-aggOffset,\
                node[0,1]-aggDiameter/2-maxAggD/2-aggOffset,node[0,2]-\
                aggDiameter/2-maxAggD/2-aggOffset]))
            binMax = np.array(([node[0,0]+aggDiameter/2+maxAggD/2+aggOffset,\
                node[0,1]+aggDiameter/2+maxAggD/2+aggOffset,node[0,2]+\
                aggDiameter/2+maxAggD/2+aggOffset]))


            # Check if particle overlapping any existing particles or bad nodes
            overlap = self.overlapCheck(node,aggDiameter,trianglePoints,binMin,\
                binMax,minAggD,maxEdgeLength,aggOffset,aggDiameterList)

            if overlap[0] == False:

                # Temporarily set this and other instances to check regardless if critical
                if overlap[1] == True or overlap[1] == False:

                    # Check if particle is inside the mesh if critically close          
                    inside = self.insideCheck(vertices,\
                    tets,node,aggDiameter,binMin,binMax,coord1,coord2,coord3,\
                    coord4,maxC)

                else:

                    inside = True

                # Indicate placed particle and break While Loop
                if inside == True and overlap[0] == False:
                    ParticleGen.node = node
                    break

        if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
            print("%s Remaining. Attempts required: %s" % \
                (len(aggDiameterList)-x-1, int(i/3)))

        return newMaxIter

    def generateParticleMPI(self,trianglePoints,maxAggNum,minC,maxC,\
        vertices,tets,coord1,coord2,coord3,coord4,newMaxIter,maxIter,minAggD,maxAggD,\
        aggOffset,verbose,aggDiameterList,maxEdgeLength,aggDiameter):
        
        # Generate random numbers to use in generation
        randomN = np.random.rand(newMaxIter*3)
        i=0
        ntet = len(tets)
        # Generate random nodal location
        while True:
            i=i+3

            if i/3 >= newMaxIter:
                i = 3
                newMaxIter = newMaxIter*2
                randomN = np.random.rand(newMaxIter*3)

            if newMaxIter >= maxIter:
                print("This particle has exceeeded the %r specified maximum iterations allowed." % (maxIter))
                print('Now exitting...')
                exit()

            # Random point selection in random tet prism container
            tetVerts = np.vstack((vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),0]-1),:],\
                vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),1]-1),:],\
                vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),2]-1),:],\
                vertices[int(tets[int(int(np.around(randomN[i]*ntet))-1),3]-1),:]))

            tetMin = np.amin(tetVerts, axis=0)
            tetMax = np.amax(tetVerts, axis=0)

            node = np.array([randomN[i]*(tetMax[0]-tetMin[0])+tetMin[0],\
                randomN[i+1]*(tetMax[1]-tetMin[1])+tetMin[1],randomN[i+2]\
                *(tetMax[2]-tetMin[2])+tetMin[2]]).T
            node = node[np.newaxis,:]           


            # Obtain extents for floating bin
            binMin = np.array(([node[0,0]-aggDiameter/2-maxAggD/2-aggOffset,\
                node[0,1]-aggDiameter/2-maxAggD/2-aggOffset,node[0,2]-\
                aggDiameter/2-maxAggD/2-aggOffset]))
            binMax = np.array(([node[0,0]+aggDiameter/2+maxAggD/2+aggOffset,\
                node[0,1]+aggDiameter/2+maxAggD/2+aggOffset,node[0,2]+\
                aggDiameter/2+maxAggD/2+aggOffset]))


            # Check if particle overlapping any existing particles or bad nodes
            overlap = self.overlapCheck(node,aggDiameter,trianglePoints,binMin,\
                binMax,minAggD,maxEdgeLength,aggOffset,aggDiameterList)

            if overlap[0] == False:

                
                if overlap[1] == True or overlap[1] == False:


                    # Check if particle is inside the mesh if critically close          
                    inside = self.insideCheck(vertices,\
                    tets,node,aggDiameter,binMin,binMax,coord1,coord2,coord3,\
                    coord4,maxC)

                else:

                    inside = True

                # Indicate placed particle and break While Loop
                if inside == True and overlap[0] == False:
                    ParticleGen.node = node
                    break


        return np.append(node[0,:],[aggDiameter,newMaxIter,int(i/3)])


    def generateSubParticle(self,x,trianglePoints,aggDiameter,maxAggNum,minC,\
        maxC,vertices,tets,coord1,coord2,coord3,coord4,maxIter,minAggD,maxAggD,\
        aggOffset,verbose,aggDiameterList,voxels,X_Size,Y_Size,Z_Size,\
        resolution,aggGrainsDiameterList, binderDiameterList, itzDiameterList,\
        material,maxEdgeLength):
        
        # Generate random numbers to use in generation
        randomN = np.random.rand(maxIter*3)
        i=0
        nVoxels = len(voxels)

        # Generate random nodal location
        while True:
            i=i+3

            if i >= len(randomN):
                print("This particle has exceeeded the %r specified maximum iterations allowed." % (maxIter))
                print('Now exitting...')
                exit()

            # Selection of random voxel
            randVoxel = round(randomN[i]*nVoxels)

            # Find voxel coordinates of random voxel
            xVoxel = np.floor(voxels[randVoxel-1]/(Z_Size*Y_Size))+1
            yVoxel = np.floor((voxels[randVoxel-1]-(xVoxel-1)*(Z_Size*Y_Size))/Z_Size)+1
            zVoxel = voxels[randVoxel-1]-(xVoxel-1)*(Z_Size*Y_Size)-(yVoxel-1)*Z_Size

            # Min and max coordinates of voxel (offset to align with geometry)
            voxMin = np.array(((xVoxel-1)*resolution-resolution/2,(yVoxel-1)*\
                resolution-resolution/2,(zVoxel-1)*resolution-resolution/2))+minC
            voxMax = np.array(((xVoxel-1)*resolution+resolution/2,(yVoxel-1)*\
                resolution+resolution/2,(zVoxel-1)*resolution+resolution/2))+minC

            # Check voxel max is within max of geometry
            if (voxMax[0]<maxC[0] and voxMax[1]<maxC[1] and voxMax[2]<maxC[2]).all(): 

                # Select random point in voxel
                node = np.array([randomN[i]*(voxMax[0]-voxMin[0])+voxMin[0],\
                    randomN[i+1]*(voxMax[1]-voxMin[1])+voxMin[1],randomN[i+2]\
                    *(voxMax[2]-voxMin[2])+voxMin[2]]).T
                node = node[np.newaxis,:]           

                # Obtain extents for floating bin
                binMin = np.array(([node[0,0]-aggDiameter/2-maxAggD/2-\
                    aggOffset,node[0,1]-aggDiameter/2-maxAggD/2-aggOffset,\
                    node[0,2]-aggDiameter/2-maxAggD/2-aggOffset]))
                binMax = np.array(([node[0,0]+aggDiameter/2+maxAggD/2+\
                    aggOffset,node[0,1]+aggDiameter/2+maxAggD/2+aggOffset,\
                    node[0,2]+aggDiameter/2+maxAggD/2+aggOffset]))


                # Check if particle overlapping any existing particles or bad nodes
                overlap = self.overlapCheck(node,aggDiameter,trianglePoints,\
                    binMin,binMax,minAggD,maxEdgeLength,aggOffset,np.concatenate((\
                    aggGrainsDiameterList, itzDiameterList, binderDiameterList)))


                if overlap[0] == False:

                    
                    if overlap[1] == True or overlap[1] == False:


                        # Check if particle is inside the mesh if critically close          
                        inside = self.insideCheck(vertices,\
                        tets,node,aggDiameter,binMin,binMax,coord1,coord2,coord3,\
                        coord4,maxC)

                    else:

                        inside = True

                    # Indicate placed particle and break While Loop
                    if inside == True and overlap[0] == False:
                        ParticleGen.node = node
                        break


        if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
            print("%s %s Grains Remaining. Attempts required: %s" % (
                len(aggDiameterList)-x, material, int(i/3)))

    def generateSubParticleCement(self,x,trianglePoints,aggDiameter,maxAggNum,minC,\
        maxC,vertices,tets,coord1,coord2,coord3,coord4,maxIter,minAggD,maxAggD,\
        aggOffset,verbose,aggDiameterList,voxels,CementStructure_X_Size,CementStructure_Y_Size,\
        CementStructure_Z_Size,cementStructureResolution,PoresDiameterList,ClinkerDiameterList,\
        CHDiameterList,CSH_LDDiameterList,CSH_HDDiameterList,material,maxEdgeLength):
        
        # Generate random numbers to use in generation
        randomN = np.random.rand(maxIter*3)
        i=0
        nVoxels = len(voxels)
        #print(voxels)

        # Generate random nodal location
        while True:
            i=i+3

            if i >= len(randomN):
                print("This particle has exceeeded the %r specified maximum iterations allowed." % (maxIter))
                print('Now exitting...')
                exit()

            # Selection of random voxel
            randVoxel = np.int(round(randomN[i]*nVoxels))

            # Find voxel coordinates of random voxel
            xVoxel = np.floor(voxels[randVoxel-1]/(CementStructure_Z_Size*CementStructure_Y_Size))+1
            yVoxel = np.floor((voxels[randVoxel-1]-(xVoxel-1)*(CementStructure_Z_Size*CementStructure_Y_Size))/CementStructure_Z_Size)+1
            zVoxel = voxels[randVoxel-1]-(xVoxel-1)*(CementStructure_Z_Size*CementStructure_Y_Size)-(yVoxel-1)*CementStructure_Z_Size+1
            #zVoxel = voxels[randVoxel-1]-(xVoxel-1)*(CementStructure_Z_Size*CementStructure_Y_Size)-(yVoxel-1)*CementStructure_Z_Size

            # Min and max coordinates of voxel (offset to align with geometry)
            voxMin = np.array(((xVoxel-1)*cementStructureResolution,(yVoxel-1)*\
                cementStructureResolution,(zVoxel-1)*cementStructureResolution))+minC
            voxMax = np.array(((xVoxel)*cementStructureResolution,(yVoxel)*\
                cementStructureResolution,(zVoxel-1)*cementStructureResolution))+minC

            #voxMin = np.array(((xVoxel-1)*cementStructureResolution-cementStructureResolution/2,(yVoxel-1)*\
            #    cementStructureResolution-cementStructureResolution/2,(zVoxel-1)*cementStructureResolution-cementStructureResolution/2))+minC
            #voxMax = np.array(((xVoxel-1)*cementStructureResolution+cementStructureResolution/2,(yVoxel-1)*\
            #    cementStructureResolution+cementStructureResolution/2,(zVoxel-1)*cementStructureResolution+cementStructureResolution/2))+minC

            # Check voxel max is within max of geometry
            if (voxMax[0]<maxC[0] and voxMax[1]<maxC[1] and voxMax[2]<maxC[2]).all(): 

                # Select random point in voxel
                node = np.array([randomN[i]*(voxMax[0]-voxMin[0])+voxMin[0],\
                    randomN[i+1]*(voxMax[1]-voxMin[1])+voxMin[1],randomN[i+2]\
                    *(voxMax[2]-voxMin[2])+voxMin[2]]).T
                node = node[np.newaxis,:] 
                #print(node)          

                # Obtain extents for floating bin
                binMin = np.array(([node[0,0]-aggDiameter/2-maxAggD/2-\
                    aggOffset,node[0,1]-aggDiameter/2-maxAggD/2-aggOffset,\
                    node[0,2]-aggDiameter/2-maxAggD/2-aggOffset]))
                binMax = np.array(([node[0,0]+aggDiameter/2+maxAggD/2+\
                    aggOffset,node[0,1]+aggDiameter/2+maxAggD/2+aggOffset,\
                    node[0,2]+aggDiameter/2+maxAggD/2+aggOffset]))


                # Check if particle overlapping any existing particles or bad nodes
                overlap = self.overlapCheck(node,aggDiameter,trianglePoints,\
                    binMin,binMax,minAggD,maxEdgeLength,aggOffset,np.concatenate((\
                    PoresDiameterList, ClinkerDiameterList, CHDiameterList,\
                    CSH_LDDiameterList,CSH_HDDiameterList)))

                if overlap[0] == False:

                    
                    if overlap[1] == True or overlap[1] == False:


                        # Check if particle is inside the mesh if critically close          
                        inside = self.insideCheck(vertices,\
                        tets,node,aggDiameter,binMin,binMax,coord1,coord2,coord3,\
                        coord4,maxC)

                    else:

                        inside = True

                    # Indicate placed particle and break While Loop
                    if inside == True and overlap[0] == False:
                        ParticleGen.node = node
                        break

        if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
            print("%s %s Grains Remaining. Attempts required: %s" % (
                len(aggDiameterList)-x, material, int(i/3)))

    def overlapCheckMPI(self,center,aggDiameter,binMin,binMax,minAggD,aggOffset,nodes,diameters):

        # Store particle nodes that fall inside the bin
        binTestParticles = np.all([(nodes[:,0] > binMin[0]) , \
            (nodes[:,0] < binMax[0]) , (nodes[:,1] > binMin[1]) , \
            (nodes[:,1] < binMax[1]) , (nodes[:,2] > binMin[2]) , \
            (nodes[:,2] < binMax[2])],axis=0)
        existingNodes = nodes[binTestParticles,:]
        existingAggD = diameters[binTestParticles]

        # Compute distance between particles 
        if len(existingNodes>0):
            nodalDistance = np.linalg.norm(center-existingNodes, axis=1)
            aggOffsetDist = nodalDistance - aggDiameter/2 - existingAggD\
                /2 - aggOffset
        else: 
            aggOffsetDist = np.array(([1]))

        # Kill and return if overlap
        if (aggOffsetDist<0).any():
            return True

        # Otherwise return false
        return False


    def overlapCheck(self,center,aggDiameter,trianglePoints,binMin,binMax,\
        minAggD,maxEdgeLength,aggOffset,aggDiameterList):

        # Store particle nodes that fall inside the bin
        binTestParticles = np.all([(self.nodes[:,0] > binMin[0]) , \
            (self.nodes[:,0] < binMax[0]) , (self.nodes[:,1] > binMin[1]) , \
            (self.nodes[:,1] < binMax[1]) , (self.nodes[:,2] > binMin[2]) , \
            (self.nodes[:,2] < binMax[2])],axis=0)
        existingNodes = self.nodes[binTestParticles,:]
        existingAggD = aggDiameterList[binTestParticles]

        # Compute distance between particles 
        if len(existingNodes>0):
            nodalDistance = np.linalg.norm(center-existingNodes, axis=1)
            aggOffsetDist = nodalDistance - aggDiameter/2 - existingAggD\
                /2 - aggOffset
        else: 
            aggOffsetDist = np.array(([1]))

        # Kill and return if overlap
        if (aggOffsetDist<0).any():
            return True,"NA"

        # Store edge nodes that fall inside the bin
        binTestSurf = np.all([(trianglePoints[:,0] > binMin[0]) , \
            (trianglePoints[:,0] < binMax[0]) , (trianglePoints[:,1] > \
                binMin[1]) , (trianglePoints[:,1] < binMax[1]) ,\
            (trianglePoints[:,2] > binMin[2]) , (trianglePoints[:,2] < \
                binMax[2])],axis=0)     
        existingSurf = trianglePoints[binTestSurf,:]

        # Compute distance between particle and edge nodes
        if len(existingSurf>0):
            surfNodalDistance = np.linalg.norm(center-existingSurf, axis=1)
            aggSurfaceDist = surfNodalDistance - aggDiameter/2 - 1.1*minAggD/2
        else: 
            aggSurfaceDist = np.array(([1]))

        # Kill and return if overlap
        if (aggSurfaceDist<0).any():
            return True,"NA"

        # Otherwise return false and check if critically close to surface
        if len(existingSurf>0):
            if (surfNodalDistance<math.sqrt(1/3*maxEdgeLength**2+(aggDiameter/2)**2)).any():
                return False,True

        return False,False


    # Checks whether all particle voxel nodes fall completely within tets 
    # (thus inside). 
    def insideCheck(self,vertices,tets,center,aggDiameter,binMin,binMax,coord1,\
        coord2,coord3,coord4,maxC):

        # Store tet vertices that fall inside the bin
        coord1 = np.all([(coord1[:,0] > binMin[0]) , (coord1[:,0] < binMax[0]),\
            (coord1[:,1] > binMin[1]) , (coord1[:,1] < binMax[1]) ,\
            (coord1[:,2] > binMin[2]) , (coord1[:,2] < binMax[2])],axis=0)      
        coord2 = np.all([(coord2[:,0] > binMin[0]) , (coord2[:,0] < binMax[0]),\
            (coord2[:,1] > binMin[1]) , (coord2[:,1] < binMax[1]) ,\
            (coord2[:,2] > binMin[2]) , (coord2[:,2] < binMax[2])],axis=0)          
        coord3 = np.all([(coord3[:,0] > binMin[0]) , (coord3[:,0] < binMax[0]),\
            (coord3[:,1] > binMin[1]) , (coord3[:,1] < binMax[1]) ,\
            (coord3[:,2] > binMin[2]) , (coord3[:,2] < binMax[2])],axis=0)  
        coord4 = np.all([(coord4[:,0] > binMin[0]) , (coord4[:,0] < binMax[0]),\
            (coord4[:,1] > binMin[1]) , (coord4[:,1] < binMax[1]) ,\
            (coord4[:,2] > binMin[2]) , (coord4[:,2] < binMax[2])],axis=0)  

        binTets = np.any([coord1,coord2,coord3,coord4],axis=0)

        coord1 = vertices[tets.astype(int)[binTets,0]-1]
        coord2 = vertices[tets.astype(int)[binTets,1]-1]
        coord3 = vertices[tets.astype(int)[binTets,2]-1]
        coord4 = vertices[tets.astype(int)[binTets,3]-1]

        node = np.empty((8,3))

        node[0,:] = [center[0,0]+aggDiameter/2,center[0,1]+aggDiameter/2,\
            center[0,2]+aggDiameter/2]
        node[1,:] = [center[0,0]-aggDiameter/2,center[0,1]+aggDiameter/2,\
            center[0,2]+aggDiameter/2]
        node[2,:] = [center[0,0]+aggDiameter/2,center[0,1]-aggDiameter/2,\
            center[0,2]+aggDiameter/2]
        node[3,:] = [center[0,0]+aggDiameter/2,center[0,1]+aggDiameter/2,\
            center[0,2]-aggDiameter/2]
        node[4,:] = [center[0,0]-aggDiameter/2,center[0,1]-aggDiameter/2,\
            center[0,2]+aggDiameter/2]
        node[5,:] = [center[0,0]+aggDiameter/2,center[0,1]-aggDiameter/2,\
            center[0,2]-aggDiameter/2]
        node[6,:] = [center[0,0]-aggDiameter/2,center[0,1]+aggDiameter/2,\
            center[0,2]-aggDiameter/2]
        node[7,:] = [center[0,0]-aggDiameter/2,center[0,1]-aggDiameter/2,\
            center[0,2]-aggDiameter/2]

        inside = 0
        emptyOnes = np.ones(len(coord1[:,0]))

        D00 = np.rot90(np.dstack((coord1[:,0],coord1[:,1],coord1[:,2],\
            emptyOnes)), 3)
        D01 = np.rot90(np.dstack((coord2[:,0],coord2[:,1],coord2[:,2],\
            emptyOnes)), 3)
        D02 = np.rot90(np.dstack((coord3[:,0],coord3[:,1],coord3[:,2],\
            emptyOnes)), 3)
        D03 = np.rot90(np.dstack((coord4[:,0],coord4[:,1],coord4[:,2],\
            emptyOnes)), 3)

        D0 = np.linalg.det(np.hstack((D00,D01,D02,D03)))

        for i in range(8):

            D10 = np.rot90(np.dstack((emptyOnes*node[i,0],\
                emptyOnes*node[i,1],emptyOnes*node[i,2],emptyOnes)), 3)
            
            D1 = np.linalg.det(np.hstack((D10,D01,D02,D03)))
            D2 = np.linalg.det(np.hstack((D00,D10,D02,D03)))
            D3 = np.linalg.det(np.hstack((D00,D01,D10,D03)))
            D4 = np.linalg.det(np.hstack((D00,D01,D02,D10)))

            if np.logical_and(np.logical_and(np.sign(D0) == np.sign(D1),\
                np.sign(D0) == np.sign(D2)),\
                np.logical_and(np.sign(D0) == np.sign(D3),\
                np.sign(D0) == np.sign(D4))).any():
                inside = inside+1
            else:
                #if ParticleGen.badItem % 100 == 0:
                #    ParticleGen.badList = np.concatenate((ParticleGen.\
                #        badList,np.array([maxC,]*100)*2))

                #ParticleGen.badList[ParticleGen.badItem,:] = node[i,:]
                #ParticleGen.badItem = ParticleGen.badItem+1
                pass

            if inside <= i:
                return False
                break

            if inside == 8:
                return True
                break

    def CTScanFibers(self,ctData):
        
        fibersList = np.loadtxt(ctData, usecols=range(7), \
            skiprows=3)
        CTScanFiberData=[]

        for i in range(len(fibersList)):
            p1Fiber = fibersList[i,1:4]
            p2Fiber = fibersList[i,4:7]

            orienFiber = (p2Fiber-p1Fiber)/\
                np.linalg.norm(p1Fiber-p2Fiber)
            fiberLength = np.linalg.norm(p1Fiber-p2Fiber)

            Data = np.array([p1Fiber[0],p1Fiber[1],p1Fiber[2],p2Fiber[0],p2Fiber[1],p2Fiber[2],\
                orienFiber[0],orienFiber[1],orienFiber[2],fiberLength])
            CTScanFiberData.append(Data)

        CTScanFiberData = np.array(CTScanFiberData).reshape(-1,10)

        return CTScanFiberData
        
    def generateFibers(self,vertices,tets,coord1,coord2,coord3,coord4,maxIter,\
        lFiber,maxC,maxAggD,fiberOrientation,orientationStrength,triangles,\
        cutFiber):
        
        # Generate random numbers to use in generation
        randomN1 = np.random.rand(maxIter*3)
        randomN2 = np.random.rand(maxIter*3)
        i=0
        ntet = len(tets)
        # Generate random nodal location
        while True:
            i=i+3

            if i >= len(randomN1):
                print("This fiber has exceeeded the %r specified maximum iterations allowed." % (maxIter))
                print('Now exitting...')
                exit()

            # Random point selection in random tet prism container
            tetVerts = np.vstack((vertices[tets[round(randomN1[i]*\
                ntet).astype(int)-1,0].astype(int)-1,:],\
                vertices[tets[round(randomN1[i]*ntet).astype(int)-1,1]\
                .astype(int)-1,:],\
                vertices[tets[round(randomN1[i]*ntet).astype(int)-1,2]\
                .astype(int)-1,:],\
                vertices[tets[round(randomN1[i]*ntet).astype(int)-1,3]\
                .astype(int)-1,:]))

            tetMin = np.amin(tetVerts, axis=0)
            tetMax = np.amax(tetVerts, axis=0)

            p1Fiber = np.array([randomN1[i]*(tetMax[0]-tetMin[0])+tetMin[0],\
                randomN1[i+1]*(tetMax[1]-tetMin[1])+tetMin[1],randomN1[i+2]\
                *(tetMax[2]-tetMin[2])+tetMin[2]]).T        

            if fiberOrientation == []:

                # Option for Totally Random Orientation (Get spherical -> Cartesian -> Normalize)

                orienFiber1 = np.array((1,randomN2[i+1]*2*np.pi,randomN2[i+2]*np.pi))

                orienFiber2 = np.array((np.sin(orienFiber1[2])*np.cos(orienFiber1[1]),np.sin(orienFiber1[2])*np.sin(orienFiber1[1]),np.cos(orienFiber1[2])))

                orienFiber = np.array((orienFiber2[0],orienFiber2[1],orienFiber2[2]))/\
                    np.linalg.norm(np.array((orienFiber2[0],orienFiber2[1]\
                    ,orienFiber2[2])))

            else:

                # Option with Preferred Orientation

                strength = (6**(4-4*orientationStrength)-1)/200
            
                v = np.empty(2)
                j = 0

                while j < 2:
                    y = np.random.normal(0, strength, 1)
                    if y > -1 and y < 1:
                        v[j] = y
                        j = j+1

                # Normalize fiber orientation
                orienFiber1 = np.array((fiberOrientation[0],fiberOrientation[1],fiberOrientation[2]))/\
                    np.linalg.norm(np.array((fiberOrientation[0],fiberOrientation[1],fiberOrientation[2])))

                # Get spherical coordinates
                orienFiber2 = np.array((1,np.arctan2(orienFiber1[1],orienFiber1[0]),np.arccos(orienFiber1[2]/(orienFiber1[0]**2+orienFiber1[1]**2+orienFiber1[2]**2)**0.5)))

                # Perturb values
                orienFiber3 = np.array((1,orienFiber2[1]+np.pi*v[0],orienFiber2[2]+np.pi/2*v[1]))

                # Convert back to Cartesian
                orienFiber = np.array((np.sin(orienFiber3[2])*np.cos(orienFiber3[1]),np.sin(orienFiber3[2])*np.sin(orienFiber3[1]),np.cos(orienFiber3[2])))

                randSign = np.random.rand(1)
                if randSign<0.5:
                    sign = -1
                else:
                    sign = 1

                # Include opposite direction
                orienFiber = sign*orienFiber


            p2Fiber = p1Fiber+orienFiber*lFiber

            # Obtain extents for floating bin
            binMin = np.amin(np.vstack((p1Fiber,p2Fiber)), axis=0)-maxAggD-lFiber
            binMax = np.amax(np.vstack((p1Fiber,p2Fiber)), axis=0)+maxAggD+lFiber

            # Check if fiber is inside the mesh    
            inside = False     
            inside = self.insideCheckFiber(vertices,tets,p1Fiber,p2Fiber,\
                binMin,binMax,coord1,coord2,coord3,coord4,maxC)

            # Indicate placed fiber and break While Loop
            if inside == True:
                ParticleGen.p1Fiber = p1Fiber
                ParticleGen.p2Fiber = p2Fiber
                ParticleGen.orienFiber = orienFiber
                ParticleGen.fiberLength = lFiber
                break
            
            # Find point fiber intersects external surface and trim accordingly
            else:

                # Find point fiber intersects external surface and trim accordingly
                if cutFiber in ['on','On','Y','y','Yes','yes']:                  

                    # Get all surface triangle coordinates
                    triangles = triangles.astype(int)
                
                    coords0 = vertices[triangles[:,0]-1]
                    coords1 = vertices[triangles[:,1]-1]
                    coords2 = vertices[triangles[:,2]-1] 

                    averageTriangles = (coords0+coords1+coords2)/3
                    averageFiber = (p1Fiber+p2Fiber)/2

                    # Find distance to nearest surface triangle
                    distances = np.linalg.norm(averageTriangles-p2Fiber,axis=1)
                    nearest = np.where(distances == np.amin(distances))

                    # Store the plane of this triangle
                    p0 = coords0[nearest,:]
                    p1 = coords1[nearest,:]
                    p2 = coords2[nearest,:]

                    p01 = p1-p0
                    p02 = p2-p0

                    fiberVector = p2Fiber-p1Fiber
 
                    # Compute distance to cutting plane
                    t = (np.dot(np.squeeze(np.cross(p01,p02)),np.squeeze((p1Fiber-p0))))/(np.dot(np.squeeze(-fiberVector),np.squeeze(np.cross(p01,p02))))

                   # New point 2 for fiber after cutting
                    p2Fiber = p1Fiber+fiberVector*t

                   # Obtain extents for floating bin
                    binMin = np.amin(np.vstack((p1Fiber,p2Fiber)), axis=0)-maxAggD-np.linalg.norm(p1Fiber-p2Fiber)
                    binMax = np.amax(np.vstack((p1Fiber,p2Fiber)), axis=0)+maxAggD+np.linalg.norm(p1Fiber-p2Fiber)

                    # Verfiy cut fiber is inside the mesh      
                    inside = False   
                    inside = self.insideCheckFiber(vertices,tets,p1Fiber,p1Fiber+0.99999*fiberVector*t,\
                        binMin,binMax,coord1,coord2,coord3,coord4,maxC)

                    if np.logical_and(inside == True,np.linalg.norm(p1Fiber-p2Fiber)<lFiber):
                        ParticleGen.p2Fiber = p2Fiber
                        ParticleGen.p1Fiber = p1Fiber
                        ParticleGen.orienFiber = orienFiber
                        ParticleGen.fiberLength = np.linalg.norm(p1Fiber-p2Fiber)
                        break

                # If not trimming then discard fiber and try again
                else:

                    pass

    # Checks whether second fiber node falls within a tet (thus inside). 
    def insideCheckFiber(self,vertices,tets,p1Fiber,p2Fiber,binMin,binMax,\
        coord1,coord2,coord3,coord4,maxC):

        # Store tet vertices that fall inside the bin
        coord1 = np.all([(coord1[:,0] > binMin[0]) , (coord1[:,0] < binMax[0]),\
            (coord1[:,1] > binMin[1]) , (coord1[:,1] < binMax[1]) ,\
            (coord1[:,2] > binMin[2]) , (coord1[:,2] < binMax[2])],axis=0)      
        coord2 = np.all([(coord2[:,0] > binMin[0]) , (coord2[:,0] < binMax[0]),\
            (coord2[:,1] > binMin[1]) , (coord2[:,1] < binMax[1]) ,\
            (coord2[:,2] > binMin[2]) , (coord2[:,2] < binMax[2])],axis=0)          
        coord3 = np.all([(coord3[:,0] > binMin[0]) , (coord3[:,0] < binMax[0]),\
            (coord3[:,1] > binMin[1]) , (coord3[:,1] < binMax[1]) ,\
            (coord3[:,2] > binMin[2]) , (coord3[:,2] < binMax[2])],axis=0)  
        coord4 = np.all([(coord4[:,0] > binMin[0]) , (coord4[:,0] < binMax[0]),\
            (coord4[:,1] > binMin[1]) , (coord4[:,1] < binMax[1]) ,\
            (coord4[:,2] > binMin[2]) , (coord4[:,2] < binMax[2])],axis=0)  

        binTets = np.any([coord1,coord2,coord3,coord4],axis=0)

        coord1 = vertices[tets.astype(int)[binTets,0]-1]
        coord2 = vertices[tets.astype(int)[binTets,1]-1]
        coord3 = vertices[tets.astype(int)[binTets,2]-1]
        coord4 = vertices[tets.astype(int)[binTets,3]-1]

        emptyOnes = np.ones(len(coord1[:,0]))

        D00 = np.rot90(np.dstack((coord1[:,0],coord1[:,1],coord1[:,2],\
            emptyOnes)), 3)
        D01 = np.rot90(np.dstack((coord2[:,0],coord2[:,1],coord2[:,2],\
            emptyOnes)), 3)
        D02 = np.rot90(np.dstack((coord3[:,0],coord3[:,1],coord3[:,2],\
            emptyOnes)), 3)
        D03 = np.rot90(np.dstack((coord4[:,0],coord4[:,1],coord4[:,2],\
            emptyOnes)), 3)

        D0 = np.linalg.det(np.hstack((D00,D01,D02,D03)))

        D10 = np.rot90(np.dstack((emptyOnes*p1Fiber[0],\
            emptyOnes*p1Fiber[1],emptyOnes*p1Fiber[2],emptyOnes)), 3)
        
        D1 = np.linalg.det(np.hstack((D10,D01,D02,D03)))
        D2 = np.linalg.det(np.hstack((D00,D10,D02,D03)))
        D3 = np.linalg.det(np.hstack((D00,D01,D10,D03)))
        D4 = np.linalg.det(np.hstack((D00,D01,D02,D10)))

        p1 = False

        if np.logical_and(np.logical_and(np.sign(D0) == np.sign(D1),\
            np.sign(D0) == np.sign(D2)),\
            np.logical_and(np.sign(D0) == np.sign(D3),\
            np.sign(D0) == np.sign(D4))).any():
            p1 = True

        D10 = np.rot90(np.dstack((emptyOnes*p2Fiber[0],\
            emptyOnes*p2Fiber[1],emptyOnes*p2Fiber[2],emptyOnes)), 3)
        
        D1 = np.linalg.det(np.hstack((D10,D01,D02,D03)))
        D2 = np.linalg.det(np.hstack((D00,D10,D02,D03)))
        D3 = np.linalg.det(np.hstack((D00,D01,D10,D03)))
        D4 = np.linalg.det(np.hstack((D00,D01,D02,D10)))

        p2 = False

        if np.logical_and(np.logical_and(np.sign(D0) == np.sign(D1),\
            np.sign(D0) == np.sign(D2)),\
            np.logical_and(np.sign(D0) == np.sign(D3),\
            np.sign(D0) == np.sign(D4))).any():
            p2 = True

        if np.logical_and(p1 == True,p2==True):
            return True
            
    # Checks whether point falls within a tet (thus inside). 
    def insideCheckPoint(self,vertices,tets,point,binMin,binMax,\
        coord1,coord2,coord3,coord4,maxC):

        # Store tet vertices that fall inside the bin
        coord1 = np.all([(coord1[:,0] > binMin[0]) , (coord1[:,0] < binMax[0]),\
            (coord1[:,1] > binMin[1]) , (coord1[:,1] < binMax[1]) ,\
            (coord1[:,2] > binMin[2]) , (coord1[:,2] < binMax[2])],axis=0)      
        coord2 = np.all([(coord2[:,0] > binMin[0]) , (coord2[:,0] < binMax[0]),\
            (coord2[:,1] > binMin[1]) , (coord2[:,1] < binMax[1]) ,\
            (coord2[:,2] > binMin[2]) , (coord2[:,2] < binMax[2])],axis=0)          
        coord3 = np.all([(coord3[:,0] > binMin[0]) , (coord3[:,0] < binMax[0]),\
            (coord3[:,1] > binMin[1]) , (coord3[:,1] < binMax[1]) ,\
            (coord3[:,2] > binMin[2]) , (coord3[:,2] < binMax[2])],axis=0)  
        coord4 = np.all([(coord4[:,0] > binMin[0]) , (coord4[:,0] < binMax[0]),\
            (coord4[:,1] > binMin[1]) , (coord4[:,1] < binMax[1]) ,\
            (coord4[:,2] > binMin[2]) , (coord4[:,2] < binMax[2])],axis=0)  

        binTets = np.any([coord1,coord2,coord3,coord4],axis=0)

        coord1 = vertices[tets.astype(int)[binTets,0]-1]
        coord2 = vertices[tets.astype(int)[binTets,1]-1]
        coord3 = vertices[tets.astype(int)[binTets,2]-1]
        coord4 = vertices[tets.astype(int)[binTets,3]-1]

        inside = 0
        emptyOnes = np.ones(len(coord1[:,0]))

        D00 = np.rot90(np.dstack((coord1[:,0],coord1[:,1],coord1[:,2],\
            emptyOnes)), 3)
        D01 = np.rot90(np.dstack((coord2[:,0],coord2[:,1],coord2[:,2],\
            emptyOnes)), 3)
        D02 = np.rot90(np.dstack((coord3[:,0],coord3[:,1],coord3[:,2],\
            emptyOnes)), 3)
        D03 = np.rot90(np.dstack((coord4[:,0],coord4[:,1],coord4[:,2],\
            emptyOnes)), 3)

        D0 = np.linalg.det(np.hstack((D00,D01,D02,D03)))

        D10 = np.rot90(np.dstack((emptyOnes*point[0],\
            emptyOnes*point[1],emptyOnes*point[2],emptyOnes)), 3)
        
        D1 = np.linalg.det(np.hstack((D10,D01,D02,D03)))
        D2 = np.linalg.det(np.hstack((D00,D10,D02,D03)))
        D3 = np.linalg.det(np.hstack((D00,D01,D10,D03)))
        D4 = np.linalg.det(np.hstack((D00,D01,D02,D10)))

        p = False

        if np.logical_and(np.logical_and(np.sign(D0) == np.sign(D1),\
            np.sign(D0) == np.sign(D2)),\
            np.logical_and(np.sign(D0) == np.sign(D3),\
            np.sign(D0) == np.sign(D4))).any():
            p = True

        if p == True:
            return True

    def tetrahedralization(self,vertices2D,triangles2D,geoName,geoFile,verbose):
        
        # Prepare file of internal nodes and external nodes/facets for Tetgen
        nodeRange = np.arange(len(self.nodes))+1
        nodeList = np.vstack((nodeRange,self.nodes.T)).T
        
        with open(Path('temp/' + geoName + '2D.a.node'),"w") as f:                                       
            f.write(str(len(self.nodes)) + ' 3 0 0\n ')                                   
            f.write("\n ".join(" ".join(map(str, x)) for x in nodeList))

        print('Starting Tetgen for tetrahedralization construction.')
        
        if verbose in ['O', 'o', 'On', 'on', 'Y', 'y', 'Yes', 'yes']:
            tetgenCommand = str(Path('lib/tetgen/tetgen')) + ' -pYiO0/1S0k ' \
                + str(Path('temp/' + geoName + '2D.mesh'))
        else:
            tetgenCommand = str(Path('lib/tetgen/tetgen')) + ' -pYiO0/1S0kQ ' \
                + str(Path('temp/' + geoName + '2D.mesh'))

        os.system(tetgenCommand)

        try:
            os.rename(Path('temp/' + geoName + '2D.1.vtk'),Path('temp/' + geoName \
                + '-para-mesh.000.vtk'))
        except:
            print("Tetgen may have failed during tetrahedralization. Tetgen will be started again.")
            print("If this issue persists, you may need to use another geometry or particle distribution.")
            tetGen = self.tetrahedralization(vertices2D,triangles2D,geoName,geoFile,verbose)
        os.remove(Path('temp/' + geoName + '.mesh'))
        os.remove(Path('temp/' + geoName + '2D.mesh'))
        try:
            os.remove(Path('temp/' + geoName + '2D.stl'))
        except:
            pass
        try:
            os.remove(Path('temp/' + geoName + '2D.1.edge'))
        except:
            pass
        try:
            os.remove(Path('temp/' + geoName + '2D.1.face'))
        except:
            pass
        try:    
            os.remove(Path('temp/' + geoName + '2D.a.node'))
        except:
            pass
        os.rename(Path('temp/' + geoName + '2D.1.ele'),Path('temp/' + geoName \
            + '.ele'))
        os.rename(Path('temp/' + geoName + '2D.1.node'),Path('temp/' + geoName \
            + '.node')) 
        os.replace(Path('temp/' + geoName + '-para-mesh.000.vtk'), \
            Path('meshes/' + geoName + '/' + geoName + '-para-mesh.000.vtk'))
        os.replace(Path('temp/' + geoName + '.ele'), Path('meshes/' + geoName \
            + '/' + geoName + '.ele'))
        os.replace(Path('temp/' + geoName + '.node'), Path('meshes/' + geoName \
            + '/' + geoName + '.node'))
        try:
            shutil.copyfile(Path('geometry/' + geoFile),Path( 'meshes/' +\
                geoName + '/' + geoFile))
        except:
            pass
        try:
            shutil.copyfile(Path('temp/' + geoName + '.geo'),Path( 'meshes/' +\
                geoName + '/' + geoName + '.geo'))
        except:
            pass
        try:    
            os.remove(Path('temp/' + geoName + '.geo'))
        except:
            pass


    def particleData(self,allNodes,allTets,aggDiameter,minAggD,geoName):

        # Create diameters list (including zero edge particle diameters)
        allDiameters = np.concatenate((np.array([0.0,]*\
            int(len(allNodes)-len(aggDiameter))),aggDiameter))

        # Generate data file for particles
        with open(Path('meshes/' + geoName + '/' + geoName + \
            '-data-particles.dat'),"w") as f:                                       
            f.write('# Particle Data Generated with LDPM Mesh Generation Tool\n')
            f.write('# [n x y z d]\n')
            f.write('\n')            
            f.write('# Number of Nodes: ' + str(len(allDiameters)) + '\n')    
            f.write('# Number of Aggregates: ' + str(len(aggDiameter)) + '\n')    
            for x in range(0,len(allDiameters)):
                f.write(str(x+1) + ' ' + str(allNodes[x,0]) + ' ' + str(allNodes[x,1]) \
                	+ ' ' + str(allNodes[x,2]) + ' ' + str(allDiameters[x]) + '\n')
            f.write('\n')  

        
    def vtkParticles(self,center,aggDiameter,materialList,geoName):

        # Generate VTK file for visualizing particles
        with open(Path('meshes/' + geoName + '/' + geoName + \
            '-para-particles.000.vtk'),"w") as f:                                       
            f.write('# vtk DataFile Version 2.0\n')
            f.write('Unstructured grid legacy vtk file with point scalar data\n')            
            f.write('ASCII\n')    
            f.write('\n')  
            f.write('DATASET UNSTRUCTURED_GRID\n')        
            f.write('POINTS ' + str(len(center)) + ' double \n')  
            f.write("\n".join(" ".join(map(str, x)) for x in center))
            f.write('\n\n')  
            f.write('POINT_DATA ' + str(len(center)) + '\n')
            f.write('SCALARS Diameter double\n')
            f.write('LOOKUP_TABLE default\n')
            for x in aggDiameter:
                f.write("%s\n" % x)
            f.write('\n')  
            f.write('SCALARS Material double\n')
            f.write('LOOKUP_TABLE default\n')
            for x in materialList:
                f.write("%s\n" % x)


    def vtkFibers(self,p1Fiber,p2Fiber,dFiber,lFiber,orienFibers,geoName):

        fiberCenter = (p2Fiber-p1Fiber)/2+p1Fiber
        vectorValues = np.array(([dFiber,dFiber]))
        fiberVector = np.concatenate((lFiber,np.array([vectorValues,]*len(fiberCenter))),axis=1)

        # Generate VTK file for visualizing particles
        with open(Path('meshes/' + geoName + '/' + geoName + \
            '-para-fibers.000.vtk'),"w") as f:                                       
            f.write('# vtk DataFile Version 2.0\n')
            f.write('Unstructured grid legacy vtk file with point vector data\n')            
            f.write('ASCII\n')    
            f.write('\n')  
            f.write('DATASET UNSTRUCTURED_GRID\n')        
            f.write('POINTS ' + str(len(fiberCenter)) + ' double \n')  
            f.write("\n".join(" ".join(map(str, x)) for x in p1Fiber))
            f.write('\n\n')  
            f.write('POINT_DATA ' + str(len(fiberCenter)) + '\n')
            f.write('VECTORS Size float\n')
            f.write("\n".join(" ".join(map(str, x)) for x in fiberVector))
            f.write('\n\n')  
            f.write('VECTORS Orientation float\n')
            f.write("\n".join(" ".join(map(str, x)) for x in orienFibers))

    def vtkNonIntersectedFibers(self,p1Fiber,p2Fiber,dFiber,lFiber,orienFibers,geoName,IntersectedFiber):
        Nonp1Fiber = []
        NonorienFibers = []
        fiberVector = []

        for i in range(0,len(p1Fiber)):

            if i not in IntersectedFiber:
                lenght = lFiber[i,0]
                Vector = np.array([lenght,dFiber,dFiber])
                fiberVector.append(Vector)
                Nonp1Fiber.append(p1Fiber[i])
                NonorienFibers.append(orienFibers[i])

        Nonp1Fiber=np.array(Nonp1Fiber).reshape(-1,3)
        NonorienFibers=np.array(NonorienFibers).reshape(-1,3)
        fiberVector=np.array(fiberVector).reshape(-1,3)

        # Generate VTK file for visualizing particles
        with open(Path('meshes/' + geoName + '/' + geoName + \
            '-para-nonintersectedfibers.000.vtk'),"w") as f:                                       
            f.write('# vtk DataFile Version 2.0\n')
            f.write('Unstructured grid legacy vtk file with point vector data\n')            
            f.write('ASCII\n')    
            f.write('\n')  
            f.write('DATASET UNSTRUCTURED_GRID\n')        
            f.write('POINTS ' + str(len(Nonp1Fiber)) + ' double \n')  
            f.write("\n".join(" ".join(map(str, x)) for x in Nonp1Fiber))
            f.write('\n\n')  
            f.write('POINT_DATA ' + str(len(Nonp1Fiber)) + '\n')
            f.write('VECTORS Size float\n')
            f.write("\n".join(" ".join(map(str, x)) for x in fiberVector))
            f.write('\n\n')  
            f.write('VECTORS Orientation float\n')
            f.write("\n".join(" ".join(map(str, x)) for x in NonorienFibers))

    def projectedFacetVTKFile(self,geoName,tetFacets,dataList,projectedFacet):

        facetsPoints = projectedFacet.reshape(-1,3)
        cells = (np.arange(0,round(len(facetsPoints))).\
            reshape(-1,3)).astype(int)
        cell_types = np.array([5,]*round(len(facetsPoints)/3))

        with open(Path('meshes/' + geoName + '/' + geoName + \
            '-para-projectedFacet.000.vtk'),"w") as f:                                                                          
            f.write('# vtk DataFile Version 2.0\n')
            f.write('Unstructured Grid\n')            
            f.write('ASCII\n')    
            f.write('DATASET UNSTRUCTURED_GRID\n')        
            f.write('POINTS ' + str(len(facetsPoints)) + ' double \n')  
            f.write("\n".join(" ".join(map(str, x)) for x in facetsPoints))
            f.write('\n\n')  
            f.write('CELLS ' + str(round(len(facetsPoints)/3)) + ' ' \
                + str(round(len(facetsPoints)/3*4)) +'\n3 ')
            f.write("\n3 ".join(" ".join(map(str, x)) for x in cells))
            f.write('\n\n')  
            f.write('CELL_TYPES ' + str(round(len(facetsPoints)/3)) +'\n')
            for x in cell_types:
                f.write("%s\n" % x)

    def vtkRebar(self,p1Rebar,p2Rebar,dRebar,geoName,x):

        rebarCenter = (p2Rebar-p1Rebar)/2+p1Rebar
        orienRebar = (p2Rebar-p1Rebar)
        rebarVector = np.vstack((np.linalg.norm(p2Rebar-p1Rebar, axis=1),[dRebar,]*len(rebarCenter),[dRebar,]*len(rebarCenter))).transpose()

        # Generate VTK file for visualizing particles
        with open(Path('meshes/' + geoName + '/' + geoName + \
            '-para-rebar.' + str(x+1) + '.000.vtk'),"w") as f:                                       
            f.write('# vtk DataFile Version 2.0\n')
            f.write('Unstructured grid legacy vtk file with point vector data\n')            
            f.write('ASCII\n')    
            f.write('\n')  
            f.write('DATASET UNSTRUCTURED_GRID\n')        
            f.write('POINTS ' + str(len(rebarCenter)) + ' double \n')  
            f.write("\n".join(" ".join(map(str, x)) for x in p1Rebar))
            f.write('\n\n')  
            f.write('POINT_DATA ' + str(len(rebarCenter)) + '\n')
            f.write('VECTORS Size float\n')
            for x in range(len(rebarVector)):
                f.write(str(rebarVector[x,0]) + ' ' + str(rebarVector[x,1]) + ' ' + str(rebarVector[x,2]) + ' ' + '\n')
            f.write('\n\n')  
            f.write('VECTORS Orientation float\n')
            f.write("\n".join(" ".join(map(str, x)) for x in orienRebar))


    def readTetgen(self, nodeFile, tetFile):                                       

        allNodes = np.loadtxt(nodeFile, usecols=(1,2,3), \
            skiprows=1)                                   

        allTets = np.loadtxt(tetFile, usecols=(1,2,3,4), \
            skiprows=1)   

        return allNodes, allTets
 

    def tesselation(self,allNodes,allTets,aggDiameter,minAggD,geoName):
        
        # Create diameters list (including fictitious edge particle diameters)
        allDiameters = np.concatenate((np.array([1.1*minAggD,]*\
            int(len(allNodes)-len(aggDiameter))),aggDiameter))

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

        # Definition of Face Points faceNPoint[Coordinates]

        # Face 0 (Nodes: 1,2,3) 
        faceNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,15:18], axis=1),\
            np.linalg.norm(allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,12:15], axis=1),\
            np.linalg.norm(allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,9:12], axis=1))).T
        
        faceOffsetDistance = np.vstack(((faceNodalDistance[:,0] - allDiameters[(allTets[:,1]-1).astype(int)]/2)/2,\
            (faceNodalDistance[:,1] - allDiameters[(allTets[:,2]-1).astype(int)]/2)/2,\
            (faceNodalDistance[:,2] - allDiameters[(allTets[:,3]-1).astype(int)]/2)/2)).T

        facePoints = np.hstack(((allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,15:18])/\
            np.array([faceNodalDistance[:,0],]*3).T\
            *np.array([faceOffsetDistance[:,0],]*3).T + edgePoints[:,15:18],\
            (allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,12:15])/\
            np.array([faceNodalDistance[:,1],]*3).T\
            *np.array([faceOffsetDistance[:,1],]*3).T + edgePoints[:,12:15],\
            (allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,9:12])/\
            np.array([faceNodalDistance[:,2],]*3).T\
            *np.array([faceOffsetDistance[:,2],]*3).T + edgePoints[:,9:12]))

        face1x = (facePoints[:,0]+facePoints[:,3]+facePoints[:,6])/3
        face1y = (facePoints[:,1]+facePoints[:,4]+facePoints[:,7])/3
        face1z = (facePoints[:,2]+facePoints[:,5]+facePoints[:,8])/3

        face0Point = np.vstack((face1x,face1y,face1z)).T

        # Face 1 (Nodes: 0,2,3)
        faceNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,3:6], axis=1),\
            np.linalg.norm(allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,6:9], axis=1),\
            np.linalg.norm(allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,15:18], axis=1))).T
        
        faceOffsetDistance = np.vstack(((faceNodalDistance[:,0] - allDiameters[(allTets[:,3]-1).astype(int)]/2)/2,\
            (faceNodalDistance[:,1] - allDiameters[(allTets[:,2]-1).astype(int)]/2)/2,\
            (faceNodalDistance[:,2] - allDiameters[(allTets[:,0]-1).astype(int)]/2)/2)).T

        facePoints = np.hstack(((allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,3:6])/\
            np.array([faceNodalDistance[:,0],]*3).T\
            *np.array([faceOffsetDistance[:,0],]*3).T + edgePoints[:,3:6],\
            (allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,6:9])/\
            np.array([faceNodalDistance[:,1],]*3).T\
            *np.array([faceOffsetDistance[:,1],]*3).T + edgePoints[:,6:9],\
            (allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,15:18])/\
            np.array([faceNodalDistance[:,2],]*3).T\
            *np.array([faceOffsetDistance[:,2],]*3).T + edgePoints[:,15:18]))

        face1x = (facePoints[:,0]+facePoints[:,3]+facePoints[:,6])/3
        face1y = (facePoints[:,1]+facePoints[:,4]+facePoints[:,7])/3
        face1z = (facePoints[:,2]+facePoints[:,5]+facePoints[:,8])/3

        face1Point = np.vstack((face1x,face1y,face1z)).T

        # Face 2 (Nodes: 0,1,3)
        faceNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,12:15], axis=1),\
            np.linalg.norm(allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,6:9], axis=1),\
            np.linalg.norm(allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,0:3], axis=1))).T
        
        faceOffsetDistance = np.vstack(((faceNodalDistance[:,0] - allDiameters[(allTets[:,0]-1).astype(int)]/2)/2,\
            (faceNodalDistance[:,1] - allDiameters[(allTets[:,1]-1).astype(int)]/2)/2,\
            (faceNodalDistance[:,2] - allDiameters[(allTets[:,3]-1).astype(int)]/2)/2)).T

        facePoints = np.hstack(((allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,12:15])/\
            np.array([faceNodalDistance[:,0],]*3).T\
            *np.array([faceOffsetDistance[:,0],]*3).T + edgePoints[:,12:15],\
            (allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,6:9])/\
            np.array([faceNodalDistance[:,1],]*3).T\
            *np.array([faceOffsetDistance[:,1],]*3).T + edgePoints[:,6:9],\
            (allNodes[(allTets[:,3]-1).astype(int),:]-edgePoints[:,0:3])/\
            np.array([faceNodalDistance[:,2],]*3).T\
            *np.array([faceOffsetDistance[:,2],]*3).T + edgePoints[:,0:3]))

        face1x = (facePoints[:,0]+facePoints[:,3]+facePoints[:,6])/3
        face1y = (facePoints[:,1]+facePoints[:,4]+facePoints[:,7])/3
        face1z = (facePoints[:,2]+facePoints[:,5]+facePoints[:,8])/3

        face2Point = np.vstack((face1x,face1y,face1z)).T

        # Face 3 (Nodes: 0,1,2)
        faceNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,0:3], axis=1),\
            np.linalg.norm(allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,3:6], axis=1),\
            np.linalg.norm(allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,9:12], axis=1))).T
        
        faceOffsetDistance = np.vstack(((faceNodalDistance[:,0] - allDiameters[(allTets[:,2]-1).astype(int)]/2)/2,\
            (faceNodalDistance[:,1] - allDiameters[(allTets[:,1]-1).astype(int)]/2)/2,\
            (faceNodalDistance[:,2] - allDiameters[(allTets[:,0]-1).astype(int)]/2)/2)).T

        facePoints = np.hstack(((allNodes[(allTets[:,2]-1).astype(int),:]-edgePoints[:,0:3])/\
            np.array([faceNodalDistance[:,0],]*3).T\
            *np.array([faceOffsetDistance[:,0],]*3).T + edgePoints[:,0:3],\
            (allNodes[(allTets[:,1]-1).astype(int),:]-edgePoints[:,3:6])/\
            np.array([faceNodalDistance[:,1],]*3).T\
            *np.array([faceOffsetDistance[:,1],]*3).T + edgePoints[:,3:6],\
            (allNodes[(allTets[:,0]-1).astype(int),:]-edgePoints[:,9:12])/\
            np.array([faceNodalDistance[:,2],]*3).T\
            *np.array([faceOffsetDistance[:,2],]*3).T + edgePoints[:,9:12]))

        face1x = (facePoints[:,0]+facePoints[:,3]+facePoints[:,6])/3
        face1y = (facePoints[:,1]+facePoints[:,4]+facePoints[:,7])/3
        face1z = (facePoints[:,2]+facePoints[:,5]+facePoints[:,8])/3

        face3Point = np.vstack((face1x,face1y,face1z)).T

        # Definition of Tet-Points [Coordinates]
        tetNodalDistance = np.vstack((np.linalg.norm(allNodes[(allTets[:,0]-1).astype(int),:]-face0Point, axis=1),\
            np.linalg.norm(allNodes[(allTets[:,1]-1).astype(int),:]-face1Point, axis=1),\
            np.linalg.norm(allNodes[(allTets[:,2]-1).astype(int),:]-face2Point, axis=1),\
            np.linalg.norm(allNodes[(allTets[:,3]-1).astype(int),:]-face3Point, axis=1))).T     

        tetOffsetDistance = np.vstack(((tetNodalDistance[:,0] - allDiameters[(allTets[:,0]-1).astype(int)]/2)/2,\
            (tetNodalDistance[:,1] - allDiameters[(allTets[:,1]-1).astype(int)]/2)/2,\
            (tetNodalDistance[:,2] - allDiameters[(allTets[:,2]-1).astype(int)]/2)/2,\
            (tetNodalDistance[:,3] - allDiameters[(allTets[:,3]-1).astype(int)]/2)/2)).T

        tetPoints = np.hstack(((allNodes[(allTets[:,0]-1).astype(int),:]-face0Point)/\
            np.array([tetNodalDistance[:,0],]*3).T*\
            np.array([tetOffsetDistance[:,0],]*3).T + face0Point,\
            (allNodes[(allTets[:,1]-1).astype(int),:]-face1Point)/\
            np.array([tetNodalDistance[:,1],]*3).T*\
            np.array([tetOffsetDistance[:,1],]*3).T + face1Point,\
            (allNodes[(allTets[:,2]-1).astype(int),:]-face2Point)/\
            np.array([tetNodalDistance[:,2],]*3).T*\
            np.array([tetOffsetDistance[:,2],]*3).T + face2Point,\
            (allNodes[(allTets[:,3]-1).astype(int),:]-face3Point)/\
            np.array([tetNodalDistance[:,3],]*3).T*\
            np.array([tetOffsetDistance[:,3],]*3).T + face3Point))

        tetx = (tetPoints[:,0]+tetPoints[:,3]+tetPoints[:,6]+tetPoints[:,9])/4
        tety = (tetPoints[:,1]+tetPoints[:,4]+tetPoints[:,7]+tetPoints[:,10])/4
        tetz = (tetPoints[:,2]+tetPoints[:,5]+tetPoints[:,8]+tetPoints[:,11])/4

        tetPoints = np.vstack((tetx,tety,tetz)).T

        # Coordinates for each facet in each tet
        facet1 = np.concatenate(([tetPoints,face0Point,edgePoints[:,9:12]]),axis=1)
        facet2 = np.concatenate(([tetPoints,face0Point,edgePoints[:,12:15]]),axis=1)
        facet3 = np.concatenate(([tetPoints,face0Point,edgePoints[:,15:18]]),axis=1)
        facet4 = np.concatenate(([tetPoints,face1Point,edgePoints[:,3:6]]),axis=1)
        facet5 = np.concatenate(([tetPoints,face1Point,edgePoints[:,6:9]]),axis=1)
        facet6 = np.concatenate(([tetPoints,face1Point,edgePoints[:,15:18]]),axis=1)
        facet7 = np.concatenate(([tetPoints,face2Point,edgePoints[:,0:3]]),axis=1)
        facet8 = np.concatenate(([tetPoints,face2Point,edgePoints[:,6:9]]),axis=1)
        facet9 = np.concatenate(([tetPoints,face2Point,edgePoints[:,12:15]]),axis=1)
        facet10 = np.concatenate(([tetPoints,face3Point,edgePoints[:,0:3]]),axis=1)
        facet11 = np.concatenate(([tetPoints,face3Point,edgePoints[:,3:6]]),axis=1)
        facet12 = np.concatenate(([tetPoints,face3Point,edgePoints[:,9:12]]),axis=1)

        # Facets for P0 in General Tet [facet[Point1[x,y,z],Point2[x,y,z],Point3[x,y,z]],....,facetN]
        facetsP0 = np.concatenate(([facet4,facet5,facet7,facet8,facet10,facet11]),axis=1)

        # Facets for P1 in General Tet [facet[Point1[x,y,z],Point2[x,y,z],Point3[x,y,z]],....,facetN]
        facetsP1 = np.concatenate(([facet1,facet2,facet7,facet9,facet10,facet12]),axis=1)

        # Facets for P2 in General Tet [facet[Point1[x,y,z],Point2[x,y,z],Point3[x,y,z]],....,facetN]
        facetsP2 = np.concatenate(([facet1,facet3,facet4,facet6,facet11,facet12]),axis=1)

        # Facets for P3 in General Tet [facet[Point1[x,y,z],Point2[x,y,z],Point3[x,y,z]],....,facetN]
        facetsP3 = np.concatenate(([facet2,facet3,facet5,facet6,facet8,facet9]),axis=1)

        # Combination of all facets for general tet (24)
        cellTetFacets = np.concatenate(([facetsP0,facetsP1,facetsP2,facetsP3]),\
            axis=1)

        # Combination of nonrepeating facets for general tet (12)
        tetFacets = np.concatenate(([facet1,facet2,facet3,facet4,facet5,facet6,\
            facet7,facet8,facet9,facet10,facet11,facet12]),axis=1)
        facets = tetFacets.reshape(-1,9)

        # Initialize cell facet list
        #cells = []

        #for i in range(1,len(allNodes)):

            # Find locations in tet list which contain the given vertex
        #    location = np.argwhere(allTets==i)

            # Initialize empty facet list for vertex
        #    cellFacets=np.zeros(([len(location)*6,9]))

            # Store facets for given cell
        #    for x in range(0,len(location)):
        #       cellFacets[x*6:x*6+6,:] = cellTetFacets[location[x,0],\
        #        location[x,1]*54:location[x,1]*54+54].reshape((-1,9))

            # Add cell facets to master list
        #    cells.append(cellFacets)

        threes = 3*np.array([np.ones((len(facets)))]).T

        facetCenters = np.array((facets[:,0]+facets[:,3]+facets[:,6],\
            facets[:,1]+facets[:,4]+facets[:,7],\
            facets[:,2]+facets[:,5]+facets[:,8])).T/threes

        vectorAB = np.array((tetFacets.reshape(-1,9)[:,3]-\
            tetFacets.reshape(-1,9)[:,0],tetFacets.reshape(-1,9)[:,4]-\
            tetFacets.reshape(-1,9)[:,1],tetFacets.reshape(-1,9)[:,5]-\
            tetFacets.reshape(-1,9)[:,2])).T
        vectorAC = np.array((tetFacets.reshape(-1,9)[:,6]-\
            tetFacets.reshape(-1,9)[:,0],tetFacets.reshape(-1,9)[:,7]-\
            tetFacets.reshape(-1,9)[:,1],tetFacets.reshape(-1,9)[:,8]-\
            tetFacets.reshape(-1,9)[:,2])).T

        facetAreas = 0.5 * np.linalg.norm(np.cross(vectorAB,vectorAC),axis=1)
        facetNormals = np.cross(vectorAB,vectorAC)/np.array([np.linalg.norm\
            (np.cross(vectorAB,vectorAC),axis=1),]*3).T

        # Specify tet-node connectivity for facets (i.e. facet 1 connected by 
            # node 1 and 2)
        tetn1 = [1,1,2,0,0,2,0,0,1,0,0,1]
        tetn2 = [2,3,3,2,3,3,1,3,3,1,2,2]


        return tetFacets, facetCenters, facetAreas, facetNormals, \
            tetn1, tetn2, tetPoints, allDiameters


    def edgeElementData(self,surfaceExtLength,allNodes,allTets,tetPoints,maxAggD,vertices,tets,\
            coord1,coord2,coord3,coord4,maxC):

        print('Generating edge elements...')

        # Generate list of all faces in mesh
        faces = np.sort(np.concatenate((np.vstack((allTets[:,0],allTets[:,1],allTets[:,2])).transpose(),\
            np.vstack((allTets[:,0],allTets[:,1],allTets[:,3])).transpose(),\
            np.vstack((allTets[:,0],allTets[:,2],allTets[:,3])).transpose(),\
            np.vstack((allTets[:,1],allTets[:,2],allTets[:,3])).transpose())),axis=1)

        # Generate list of unique faces in mesh
        ufaces, inverse, unique_counts = np.unique(faces,return_inverse=True,return_counts=True,axis=0)

        outerFace = np.zeros(len(ufaces))
        innerFace = np.zeros([len(ufaces),2])

        outerFaceID = np.zeros(len(ufaces))
        innerFaceID = np.zeros([len(ufaces),2])

        j = 0
        k = 0

        for i in range(0,len(ufaces)):

            location = np.where(inverse==i)[0]

            if len(location) == 1:

                outerFaceID[j] = location[0]

                if location[0] < len(faces)/4:
                    outerFace[j] = int(location[0])
                elif location[0] < 2*len(faces)/4:
                    outerFace[j] = int(location[0] - len(faces)/4)
                elif location[0] < 3*len(faces)/4:
                    outerFace[j] = int(location[0] - 2*len(faces)/4)
                else:
                    outerFace[j] = int(location[0] - 3*len(faces)/4)
                j = j+1

            else:
                for x in range(0,2):

                    innerFaceID[k,x] = location[x]

                    if location[x] < len(faces)/4:
                        innerFace[k,x] = int(location[x])
                    elif location[x] < 2*len(faces)/4:
                        innerFace[k,x] = int(location[x] - len(faces)/4)
                    elif location[x] < 3*len(faces)/4:
                        innerFace[k,x] = int(location[x] - 2*len(faces)/4)
                    else:
                        innerFace[k,x] = int(location[x] - 3*len(faces)/4)
                k = k+1

        # Trim extra rows of zeros
        outerFace = np.trim_zeros(outerFace)
        innerFace = innerFace[~np.all(innerFace == 0, axis=1)]

        outerFaceID = np.trim_zeros(outerFaceID)
        innerFaceID = innerFaceID[~np.all(innerFaceID == 0, axis=1)]

        # Vectors of sides of faces (for area calc)
        # Case 1
        c1v1 = allNodes[(faces[innerFaceID[:,0].astype(int),:][:,1]-1).astype(int),:]-allNodes[(faces[innerFaceID[:,0].astype(int),:][:,0]-1).astype(int),:]
        c1v2 = allNodes[(faces[innerFaceID[:,0].astype(int),:][:,2]-1).astype(int),:]-allNodes[(faces[innerFaceID[:,0].astype(int),:][:,0]-1).astype(int),:]
        # Case 2/3
        c2v1 = allNodes[(faces[outerFaceID.astype(int),:][:,1]-1).astype(int),:]-allNodes[(faces[outerFaceID.astype(int),:][:,0]-1).astype(int),:]
        c2v2 = allNodes[(faces[outerFaceID.astype(int),:][:,2]-1).astype(int),:]-allNodes[(faces[outerFaceID.astype(int),:][:,0]-1).astype(int),:]
        # Combined Vectors
        v1 = np.concatenate((c1v1,c2v1,c2v1))
        v2 = np.concatenate((c1v2,c2v2,c2v2))

        # Edge element points
        # Case 1: Two Tet Points
        case1Points = np.concatenate((tetPoints[(innerFace[:,0]).astype(int)],tetPoints[(innerFace[:,1]).astype(int)]),axis=1)
        # Case 2: Tet Point to Face Point
        case2Points = np.concatenate(((allNodes[(faces[outerFaceID.astype(int),:][:,0]-1).astype(int),:]+allNodes[(faces[outerFaceID.astype(int),:][:,1]-1).astype(int),:]+allNodes[(faces[outerFaceID.astype(int),:][:,2]-1).astype(int),:])/3,tetPoints[(outerFace).astype(int),:]),axis=1)
        # Case 3: Face Point to Extension
        case3Points = np.empty([len(case2Points),6])

        for x in range(len(case2Points)):
            option1 = np.concatenate((case2Points[x,0:3],surfaceExtLength*(np.cross(c2v1[x,:],c2v2[x,:]))/np.array([np.linalg.norm(np.cross(c2v1[x,:],c2v2[x,:])),]*3).T+case2Points[x,0:3]))
        
            # Obtain extents for floating bin
            binMin = option1[3:6]-maxAggD
            binMax = option1[3:6]+maxAggD

            # Check if point is inside the mesh         
            inside = self.insideCheckPoint(vertices,tets,option1[3:6],\
                binMin,binMax,coord1,coord2,coord3,coord4,maxC)

            if inside == True:
                case3Points[x,:] = np.concatenate((case2Points[x,0:3],surfaceExtLength*(np.cross(c2v2[x,:],c2v1[x,:]))/np.array([np.linalg.norm(np.cross(c2v2[x,:],c2v1[x,:])),]*3).T+case2Points[x,0:3]))

            else:
                case3Points[x,:] = option1


        # Position of internal face centers (Correction: this point should be the intersection of line T1T2 and face)
        # ref: http://geomalgorithms.com/a05-_intersect-1.html
        inCenter = np.zeros([case1Points.shape[0],3])
        for i in range(0,case1Points.shape[0]):
            normals = (np.cross(v1,v2))/np.array([np.linalg.norm(np.cross(v1,v2),axis=1),]*3).T # case 1,2 normals
            n = np.copy(normals)[i,:]
            wvector = case1Points[i,0:3] - allNodes[(faces[innerFaceID[:,0].astype(int),:][:,0]-1).astype(int),:]
            w = np.copy(wvector)[i,:]
            u = case1Points[i,3:6] - case1Points[i,0:3]
            tI = -np.dot(n,w)/np.dot(n,u)
            inCenter[i,:] = case1Points[i,0:3] + tI*u.T
       
        # inCenter = (allNodes[(faces[innerFaceID[:,0].astype(int),:][:,0]-1).astype(int),:]+allNodes[(faces[innerFaceID[:,0].astype(int),:][:,1]-1).astype(int),:]+allNodes[(faces[innerFaceID[:,0].astype(int),:][:,2]-1).astype(int),:])/3

        # Pyramidal Volume
        # Case 1, volume 1
        p1 = allNodes[(faces[innerFaceID[:,0].astype(int),:][:,0]-1).astype(int),:]
        p2 = allNodes[(faces[innerFaceID[:,0].astype(int),:][:,1]-1).astype(int),:]
        p3 = allNodes[(faces[innerFaceID[:,0].astype(int),:][:,2]-1).astype(int),:]
        p4 = case1Points[:,0:3]
        volCalc1 = np.expand_dims(np.transpose(p1-p4).T, axis=1)
        volCalc2 = np.expand_dims(np.transpose(np.cross((p2-p4),\
            (p3-p4))).T, axis=2)
        c1vol1 = np.squeeze(abs(np.matmul(volCalc1,volCalc2))/6)

        # Case 1, volume 2
        p5 = case1Points[:,3:6]
        volCalc1 = np.expand_dims(np.transpose(p1-p5).T, axis=1)
        volCalc2 = np.expand_dims(np.transpose(np.cross((p2-p5),\
            (p3-p5))).T, axis=2)
        c1vol2 = np.squeeze(abs(np.matmul(volCalc1,volCalc2))/6)

        # Case 2, total volume
        p1 = allNodes[(faces[outerFaceID.astype(int),:][:,0]-1).astype(int),:]
        p2 = allNodes[(faces[outerFaceID.astype(int),:][:,1]-1).astype(int),:]
        p3 = allNodes[(faces[outerFaceID.astype(int),:][:,2]-1).astype(int),:]
        p4 = case2Points[:,3:6]
        volCalc1 = np.expand_dims(np.transpose(p1-p4).T, axis=1)
        volCalc2 = np.expand_dims(np.transpose(np.cross((p2-p4),\
            (p3-p4))).T, axis=2)
        c2vol = np.squeeze(abs(np.matmul(volCalc1,volCalc2))/6)        

        # Case 3, total volume
        c3vol = surfaceExtLength*0.5*np.linalg.norm(np.cross(c2v1,c2v2),axis=1)

        # Edge Data Matrix
        # [x1 y1 z1 x2 y2 z2 A n1 n2 n3 L1 L2 V1 V2 c]
        edgeData = np.zeros([len(innerFace)+2*len(outerFace),15])

        edgeData[:,0:6]              = np.concatenate((case1Points,case2Points,case3Points))
        edgeData[:,6]                = 0.5*np.linalg.norm(np.cross(v1,v2),axis=1) # Area
        edgeData[:,7:10]             = (np.cross(v1,v2))/np.array([np.linalg.norm(np.cross(v1,v2),axis=1),]*3).T # case 1,2 normals
        edgeData[-case3Points.shape[0]:,7:10] = (edgeData[-case3Points.shape[0]:,0:3]-edgeData[-case3Points.shape[0]:,3:6])/np.array([np.linalg.norm(edgeData[-case3Points.shape[0]:,0:3]-edgeData[-case3Points.shape[0]:,3:6],axis=1),]*3).T # case 3 normals
        
        # Check normal is in right direction and fix if needed
        for x in range(len(edgeData)):
            edgeData[x,7:10] = edgeData[x,7:10]*((np.dot((edgeData[x,3:6]-edgeData[x,0:3]),edgeData[x,7:10]) > 0).astype(int)*2-1)
       
        edgeData[:,10]               = np.concatenate((np.linalg.norm(inCenter-case1Points[:,0:3],axis=1),np.linalg.norm(case2Points[:,0:3]-case2Points[:,3:6],axis=1),np.linalg.norm(case3Points[:,0:3]-case3Points[:,3:6],axis=1)))
       
        edgeData[0:len(inCenter),11] = np.linalg.norm(inCenter-case1Points[:,3:6],axis=1) # L2
       
        edgeData[:,12]               = np.concatenate((c1vol1,c2vol,c3vol)) # V1
       
        edgeData[0:len(inCenter),13] = c1vol2 # V2
        edgeData[:,14]               = np.concatenate((np.ones(len(case1Points)),2*np.ones(len(case2Points)),3*np.ones(len(case3Points))))

        return edgeData


    def edgeFile(self,geoName,edgeData):

        edgeFiledata = np.copy(edgeData)
        # Store Nodal Data
        nodes = np.unique(edgeFiledata[:,0:6].reshape(-1,3),axis=0)
        
        # Store Element Data
        for x in range(len(edgeFiledata)):
            
            edgeFiledata[x,0] = (np.where(~(nodes[:,0:3]-edgeFiledata[x,0:3]).any(axis=1))[0])+1
            edgeFiledata[x,1] = (np.where(~(nodes[:,0:3]-edgeFiledata[x,3:6]).any(axis=1))[0])+1     
        
        edgeFiledata = np.delete(edgeFiledata,[2,3,4,5],1)

        np.savetxt(Path('meshes/' + geoName + '/' + geoName \
            + '-data-edgeEle-geom.dat'), nodes, fmt='%.16g', delimiter=' '\
            ,header='Edge Element Data Generated with LDPM Mesh Generation Tool\n\
Number of nodes\n'+ str(len(nodes)) +
'\n\
Number of elements\n'+ str(len(edgeFiledata)) +
'\n\
[x y z]', comments='')

        with open(Path('meshes/' + geoName + '/' + geoName \
            + '-data-edgeEle-geom.dat'),'ab') as f:
            np.savetxt(f, edgeFiledata[:,:], fmt='%.16g', delimiter=' ',header='[i1 i2 A n1 n2 n3 L1 L2 V1 V2 c]'\
            , comments='')


    def edgeVTKFile(self,geoName,edgeData):

        edgeVTKdata = np.copy(edgeData) # since np.unique operation will change the original variable edgeData, use np.copy here instead
        # Store Nodal Data
        nodes = np.unique(edgeVTKdata[:,0:6].reshape(-1,3),axis=0)
        
        # Store Element Data
        for x in range(len(edgeVTKdata)):

            edgeVTKdata[x,0] = (np.where(~(nodes[:,0:3]-edgeVTKdata[x,0:3]).any(axis=1))[0])
            edgeVTKdata[x,1] = (np.where(~(nodes[:,0:3]-edgeVTKdata[x,3:6]).any(axis=1))[0])     

        edgeVTKdata = np.delete(edgeVTKdata,[2,3,4,5],1)
        
        cell_types = edgeVTKdata[:,-1].astype(int)

        cells = edgeVTKdata[:,0:2].astype(int)

        with open(Path('meshes/' + geoName + '/' + geoName + \
            '-para-edgeEle.000.vtk'),"w") as f:                                                                          
            f.write('# vtk DataFile Version 2.0\n')
            f.write('Unstructured Grid\n')            
            f.write('ASCII\n')    
            f.write('DATASET UNSTRUCTURED_GRID\n')        
            f.write('POINTS ' + str(len(nodes)) + ' double \n')  
            f.write("\n".join(" ".join(map(str, x)) for x in nodes))
            f.write('\n\n')  
            f.write('CELLS ' + str(len(edgeVTKdata)) + ' ' \
                + str(len(edgeVTKdata)*3) +'\n2 ')
            f.write("\n2 ".join(" ".join(map(str, x)) for x in cells))
            f.write('\n\n')  
            f.write('CELL_TYPES ' + str(len(edgeVTKdata)) +'\n')
            for x in cell_types:
                f.write("%s\n" % 3)


    def marsFile(self,geoName,allNodes,allTets):

        with open(Path('meshes/' + geoName + '/' + geoName + '-mesh.mrs'),\
            "w") as f:                                       
            f.write('InsertMaterial CONC LDPM-V4 {\n')
            f.write('  MixDesign {\n')
            f.write('    CementContent                  420 kg/m3\n')
            f.write('    WaterToCementRatio             0.42\n')
            f.write('    AggregateToCementRatio         4.27\n')
            f.write('    MinAggregate                   0.02 m\n')
            f.write('    MaxAggregate                   0.01 m\n')
            f.write('    FullerCoefficient              0.5\n')
            f.write('  }\n')
            f.write('  StaticParameters {\n')
            f.write('    NormalModulus                  0 Pa\n')
            f.write('    Alpha                          0.25\n')
            f.write('    TensileStrength                0 Pa\n')
            f.write('    TensileCharacteristicLength    0 m\n')
            f.write('    ShearStrengthRatio             0\n')
            f.write('    SofteningExponent              0.2\n')
            f.write('    CompressiveYieldingStrength    0 Pa\n')
            f.write('    InitialHardeningModulusRatio   0\n')
            f.write('    TransitionalStrainRatio        5\n')
            f.write('    DeviatoricStrainThresholdRatio 1\n')
            f.write('    DeviatoricDamageParameter      5\n')
            f.write('    InitialFriction                0\n')
            f.write('    AsymptoticFriction             0\n')
            f.write('    TransitionalStress             0 Pa\n')
            f.write('    DensificationRatio             1\n')
            f.write('    VolumetricDeviatoricCoupling   0\n')
            f.write('  }\n')
            f.write('}\n')
            f.write('InsertNodeList PRTC {\n')
            f.write('LengthUnits mm\n')
            f.write('Particles\n')
            f.write('ReadNodes ' + str(len(allNodes)) + '\n')  
            f.write('// i   crx    cry           crz            rad\n')
            for x in range(0,len(allNodes)):
                    f.write(str(x+1) + '    ' + str(allNodes[x,0]) + '    ' \
                        + str(allNodes[x,1]) + '    '  + str(allNodes[x,2]) \
                        + '    0\n')
            f.write('}\n') 
            f.write('ReadObjects ' + str(len(allTets)) + '\n')  
            f.write('//  i           n1          n2           n3          n4\n') 
            for x in range(0,len(allTets)):
                    f.write(str(x+1) + '    ' + str(allTets[x,0].astype(int)) \
                        + '    ' \
                        + str(allTets[x,1].astype(int)) + '    '  \
                        + str(allTets[x,2].astype(int)) \
                        + '    '  + str(allTets[x,3].astype(int)) +'\n')
            f.write('EOF')


    def abaqusFile(self,geoName,allNodes,allTets,edgeElements,edgeData,trianglePoints,triangles):

        with open(Path('meshes/' + geoName + '/' + geoName + '-mesh.inp'),\
            "w") as f:                                       

            f.write('*Heading\n')
            f.write('** Job name: ' + geoName + ' Model name: Model-' + geoName + '\n')
            f.write('** Generated by: Abaqus/CAE 2019\n')
            f.write('**\n')
            f.write('** PARTS\n')
            f.write('**\n')
            f.write('*Part, name=' + geoName + '\n')
            f.write('*Node\n')
            for x in range(0,len(allNodes)):
                    f.write(str(x+1) + ', ' + str(allNodes[x,0]) + ', ' \
                        + str(allNodes[x,1]) + ', '  + str(allNodes[x,2]) \
                        + '\n')
            f.write('*Element, type=C3D4\n')
            for x in range(0,len(allTets)):
                    f.write(str(x+1) + ', ' + str(allTets[x,0].astype(int)) \
                        + ', ' \
                        + str(allTets[x,1].astype(int)) + ', '  \
                        + str(allTets[x,2].astype(int)) \
                        + ', '  + str(allTets[x,3].astype(int)) +'\n')
            f.write('*End Part\n')
            f.write('**\n')
            f.write('**\n')
            f.write('** ASSEMBLY\n')
            f.write('**\n')
            f.write('*Assembly, name=Assembly\n')
            f.write('**\n')
            f.write('*Instance, name=' + geoName + '-1, part=' + geoName + '\n')
            f.write('*End Instance\n')
            f.write('**\n')
            f.write('*End Assembly\n')

        with open(Path('meshes/' + geoName + '/' + geoName + '-mesh2D.inp'),\
            "w") as f:                                       

            f.write('*Heading\n')
            f.write('** Job name: ' + geoName + ' Model name: Model2D-' + geoName + '\n')
            f.write('** Generated by: Abaqus/CAE 2019\n')
            f.write('**\n')
            f.write('** PARTS\n')
            f.write('**\n')
            f.write('*Part, name=' + geoName + '\n')
            f.write('*Node\n')
            for x in range(0,len(trianglePoints)):
                    f.write(str(x+1) + ', ' + str(trianglePoints[x,0]) + ', ' \
                        + str(trianglePoints[x,1]) + ', ' + str(trianglePoints[x,2]) + '\n')
            f.write('*Element, type=SFM3D3\n')
            for x in range(0,len(triangles)):
                    f.write(str(x+1) + ', ' + str(triangles[x,0].astype(int)) \
                        + ', ' \
                        + str(triangles[x,1].astype(int)) + ', '  \
                        + str(triangles[x,2].astype(int)) +'\n')
            f.write('*End Part\n')
            f.write('**\n')
            f.write('**\n')
            f.write('** ASSEMBLY\n')
            f.write('**\n')
            f.write('*Assembly, name=Assembly\n')
            f.write('**\n')
            f.write('*Instance, name=' + geoName + '-1, part=' + geoName + '\n')
            f.write('*End Instance\n')
            f.write('**\n')
            f.write('*End Assembly\n')
            
            if edgeElements in ['on','On','Y','y','Yes','yes']:
                
                edgeInpdata = np.copy(edgeData)
                
                # Store Nodal Data
                nodes = np.unique(edgeInpdata[:,0:6].reshape(-1,3),axis=0)
                
                # Store Element Data
                for x in range(len(edgeInpdata)):
        
                    edgeInpdata[x,0] = (np.where(~(nodes[:,0:3]-edgeInpdata[x,0:3]).any(axis=1))[0])+1
                    edgeInpdata[x,1] = (np.where(~(nodes[:,0:3]-edgeInpdata[x,3:6]).any(axis=1))[0])+1     
        
                edgeInpdata = np.delete(edgeInpdata,[2,3,4,5],1)
                
                nodes = np.column_stack((np.arange(0,nodes.shape[0])+1,nodes))
                nodes = nodes.astype(object)
                
                Nodeindex = np.arange(0,nodes.shape[0])+1
                nodes[:,0] = Nodeindex.astype(int)
                
                # Element Connectivity
                Connectivity = np.copy(edgeInpdata[:,0:2])
                Connectivity = np.column_stack((np.arange(0,Connectivity.shape[0])+1,Connectivity))
                Connectivity = Connectivity.astype(object)
                
                Connectivity_ghostmesh = np.copy(Connectivity)
                Connectivity_ghostmesh[:,0] = Connectivity_ghostmesh[:,0] + 10000
                
                Connectivity = Connectivity.astype(int)
                Connectivity_ghostmesh = Connectivity_ghostmesh.astype(int)
                
                # Nodes on the left
                LeftNodes = (np.where(nodes[:,3] < 0.0)[0]+1).reshape(1,-1)
                
                # Nodes on the right
                RightNodes = (np.where(nodes[:,3] >= 200.0)[0]+1).reshape(1,-1)
                
                # Nodes on the bottom
                BottomNodes = (np.where(nodes[:,2] <= 0.0)[0]+1).reshape(1,-1)
                
                # Nodes on the top
                TopNodes = (np.where(nodes[:,2] >= 100.0)[0]+1).reshape(1,-1)
                
                # Nodes on the back
                BackNodes = (np.where(nodes[:,1] <= 0.0)[0]+1).reshape(1,-1)
                
                # Nodes on the front
                FrontNodes = (np.where(nodes[:,1] >= 100.0)[0]+1).reshape(1,-1)
                
                # Monitoring Nodes 1
                MonitoringNodes1 = (np.where(abs(nodes[:,3]-30) <= 0.2)[0]+1).reshape(1,-1)
                
                # Monitoring Nodes 2
                MonitoringNodes2 = (np.where(abs(nodes[:,3]-70) <= 0.2)[0]+1).reshape(1,-1)
                
                # Monitoring Nodes 3
                MonitoringNodes3 = (np.where(abs(nodes[:,3]-120) <= 0.2)[0]+1).reshape(1,-1)
                
                # Monitoring Nodes 4
                MonitoringNodes4 = (np.where(abs(nodes[:,3]-10) <= 0.3)[0]+1).reshape(1,-1)
                
                # Monitoring Nodes 5
                MonitoringNodes5 = (np.where(abs(nodes[:,3]-50) <= 0.2)[0]+1).reshape(1,-1)
        
                
                with open(Path('meshes/' + geoName + '/' + geoName + \
                    '-data-edgeEle.inp'),"w") as f:                                                                          
                    f.write('*Heading\n') 
                    f.write('** Job name: FLE-test Model name: ' + geoName + '-data-edgeEle\n')
                    f.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')
                    f.write('**\n')
                    f.write('** PARTS\n')
                    f.write('**\n')
                    f.write('*Part, name=Part-1\n')
                    f.write('*End Part\n')
                    f.write('**\n')
                    f.write('**\n')
                    f.write('** ASSEMBLY\n')
                    f.write('**\n')
                    f.write('*Assembly, name=Assembly\n')
                    f.write('**\n')
                    f.write('*Instance, name=Part-1-1, part=Part-1\n')
                    f.write('*Node\n')
#                    f.write('*Node,Nset = AllNodes\n')
                    f.write('\n'.join(','.join(map(str, x)) for x in nodes))
                    f.write('\n')
                    f.write('*USER ELEMENT,NODES=2,TYPE=U1,PROPERTIES=20,COORDINATES=3,VARIABLES=4,UNSYMM\n')
                    f.write('12,11\n')
                    f.write('*ELEMENT,TYPE=U1,ELSET=UEL\n')
                    f.write('\n'.join(','.join(map(str, x)) for x in Connectivity))
                    f.write('\n*UEL PROPERTY, ELSET=UEL\n')
                    f.write('2.5e-6, 0.28, 5.41e-7, 2.5e-3, 1100., 4166.0, 0.05, 8.\n')
                    f.write('5000.0, 5.5, 4.0, 500000.,2.8027E-13, 3.7406E-11, 3., 2700.0\n')
                    f.write('22.85, 0.255, 1.2, 0.250\n')
                    f.write('\n')
                    f.write('*ELEMENT,TYPE=DC1D2,ELSET=VISUALIZATION\n')
                    f.write('\n'.join(','.join(map(str, x)) for x in Connectivity_ghostmesh))
                    f.write('\n')
                    f.write('\n**Section: Section-1\n')
                    f.write('*Solid Section,elset=VISUALIZATION,material=Material-1\n')
                    f.write('*End Instance\n')
                    f.write('**\n')
                    f.write('*NSET, NSET=AllNodes, instance=Part-1-1, generate\n')
                    f.write('1, ' + str(nodes[-1,0]) + ', 1\n')
                    f.write('*NSET, NSET=LeftNodes, instance=Part-1-1\n')
                    for i in range(0,round(len(LeftNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, LeftNodes[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET, NSET=RightNodes, instance=Part-1-1\n')
                    for i in range(0,round(len(RightNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, RightNodes[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET, NSET=BottomNodes, instance=Part-1-1\n')
                    for i in range(0,round(len(BottomNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, BottomNodes[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET, NSET=TopNodes, instance=Part-1-1\n')
                    for i in range(0,round(len(TopNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, TopNodes[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET, NSET=BackNodes, instance=Part-1-1\n')
                    for i in range(0,round(len(BackNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, BackNodes[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET, NSET=FrontNodes, instance=Part-1-1\n')
                    for i in range(0,round(len(FrontNodes[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, FrontNodes[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET,NSET=MonitoringNodes-1, instance=Part-1-1\n')
                    for i in range(0,round(len(MonitoringNodes1[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, MonitoringNodes1[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET,NSET=MonitoringNodes-2, instance=Part-1-1\n')
                    for i in range(0,round(len(MonitoringNodes2[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, MonitoringNodes2[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET,NSET=MonitoringNodes-3, instance=Part-1-1\n')
                    for i in range(0,round(len(MonitoringNodes3[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, MonitoringNodes3[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET,NSET=MonitoringNodes-4, instance=Part-1-1\n')
                    for i in range(0,round(len(MonitoringNodes4[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, MonitoringNodes4[0][15*i:15*(i+1)])))+'\n')
                    f.write('*NSET,NSET=MonitoringNodes-5, instance=Part-1-1\n')
                    for i in range(0,round(len(MonitoringNodes5[0])/15)): # Abaqus only accepts maximum 15 items per row
                        f.write(''.join(','.join(map(str, MonitoringNodes5[0][15*i:15*(i+1)])))+'\n')
                    f.write('*End Assembly\n')
                    f.write('**\n')
                    f.write('** MATERIALS\n')
                    f.write('**\n')
                    f.write('*Material, name=Material-1\n')
                    f.write('*Depvar\n')
                    f.write('5\n')
                    f.write('1, Phi, Phi\n')
                    f.write('2, Dh, Dh\n')
                    f.write('3, Water, Water\n')
                    f.write('4, qT, qT\n')
                    f.write('5, qH, qH\n')
                    f.write('*Density\n')
                    f.write('0.,\n')
                    f.write('*User Material, constants=1, type=THERMAL\n')
                    f.write('0.,\n')
                    f.write('**\n')
                    f.write('*Initial Conditions, type=SOLUTION\n')
                    f.write('Part-1-1.UEL, 0.3168, 0.3168\n')
                    f.write('Part-1-1.VISUALIZATION, 0.3168, 0.3168\n')
                    f.write('** PREDEFINED FIELDS\n')
                    f.write('**\n')
                    f.write('** Name: Predefined Field-1   Type: Temperature\n')
                    f.write('*Initial Conditions, type=TEMPERATURE\n')
                    f.write('AllNodes, 19.85, 0.981\n')
                    f.write('**\n')
                    f.write('**\n')
                    f.write('** STEP: Step-1\n')
                    f.write('**\n')
                    f.write('*Step, name=Step-1, nlgeom=NO, inc=9000\n')
                    f.write('*Heat Transfer, end=PERIOD\n')
                    # Initial time increment, total time period, min time increment, max time increment
                    f.write('1800, 1.62e+7 \n') 
                    f.write('**\n')
                    f.write('** BOUNDARY CONDITIONS\n')
                    f.write('**\n')
                    f.write('** Name: BC-1 Type: Temperature\n')
                    f.write('*Boundary\n')
                    f.write('LeftNodes, 11, 11, 19.85\n')
                    f.write('*Boundary\n')
                    f.write('LeftNodes, 12, 12, 0.5\n')
                    f.write('*Boundary\n')
                    f.write('RightNodes, 11, 11, 19.85\n')
                    f.write('*Boundary\n')
                    f.write('BottomNodes, 11, 11, 19.85\n')
                    f.write('*Boundary\n')
                    f.write('TopNodes, 11, 11, 19.85\n')
                    f.write('*Boundary\n')
                    f.write('BackNodes, 11, 11, 19.85\n')
                    f.write('*Boundary\n')
                    f.write('FrontNodes, 11, 11, 19.85\n')
                    f.write('**\n')
                    f.write('** OUTPUT REQUESTS\n')
                    f.write('**\n')
                    f.write('*Restart, write, frequency=10\n')
                    f.write('**\n')
                    f.write('** FIELD OUTPUT: F-Output-1\n')
                    f.write('**\n')
                    f.write('*Output, field\n')
                    f.write('*Node Output\n')
                    f.write('NT\n')
                    f.write('*Element Output, directions=YES\n')
                    f.write('SDV\n')
                    f.write('*End Step\n')


    def cellFile(self,geoName,nParticles,cells,allNodes,dataType):

        # Generate cell file for simulation

        if dataType in ['t','T','text','Text','raw','Raw']:
            with open(Path('meshes/' + geoName + '/' + geoName \
                + '-data-cell.dat'),"w") as f:                                       
                f.write('# Cells Generated with LDPM Mesh Generation Tool\n') 
                f.write('Number of LDPM Cells: ' + str(nParticles) + '\n')  
                for x in range(0,len(cells)):
                    f.write('\n\n// Cell Number\n') 
                    f.write('Cell ' + str(x+1) + '\n') 
                    f.write('// Coordinates of the center particle/node\n') 
                    f.write('crd: ' + str(allNodes[x,0]) + ' ' \
                        + str(allNodes[x,1]) + ' '  + str(allNodes[x,2]) + '\n') 
                    f.write('// Number of outside points defining cell\n') 
                    f.write('nxp: ' + str(len(cells[x])*3) + '\n') 
                    f.write('// Local coordinates of outside points follow\n') 
                    f.write("\n".join(" ".join(map(str, x)) for x in np.array(cells[x]).reshape(-1,3)))        
        else:
            np.savez(Path('meshes/' + geoName + '/' + geoName + \
                '-data-cell.dat'),cells)


    def meshFile(self,geoName,allNodes,allTets):
        
        with open(Path('meshes/' + geoName + '/' + geoName + '-data-mesh.dat'),\
            "w") as f:                                       
            f.write('# Mesh Data Generated with LDPM Mesh Generation Tool\n')
            f.write('\n') 
            f.write('Number of Vertices: ' + str(len(allNodes)) + '\n')  
            f.write("\n".join(" ".join(map(str, x)) for x in allNodes))     
            f.write('\n\n') 
            f.write('Number of Tets: ' + str(len(allTets)) + '\n')  
            f.write("\n".join(" ".join(map(str, x)) for x in \
                allTets.astype(int)))   
            f.write('\nEOF')


    def facetData(self,allNodes,allTets,tetFacets,facetCenters,\
        facetAreas,facetNormals,tetn1,tetn2,materialList,materialRule,\
        multiMaterial,cementStructure,edgeMaterialList):

        # Store facets as set of three points
        facets = tetFacets.reshape(-1,9)

        # Store particles accociated with facets
        p1 = allNodes[(allTets[:,tetn1]-1).astype(int),:].reshape(-1,3)
        p2 = allNodes[(allTets[:,tetn2]-1).astype(int),:].reshape(-1,3)


        # Projected facet normal
        pn = (p2-p1)/np.array([np.linalg.norm(p2-p1,axis=1),]*3).T
        pn = pn.reshape(-1,3)
        #print('n',pn)

        InitialNormal=[]
        coords=[]

        for x in range(0,len(allTets)):

            for y in range(0,12):
                Check = pn[12*x+y,0]*facetNormals[12*x+y,0]+\
                        pn[12*x+y,1]*facetNormals[12*x+y,1]+\
                        pn[12*x+y,2]*facetNormals[12*x+y,2]

                if Check < 0.0:
                    InitialN = -1.0*facetNormals[12*x+y,0:3]
                    c1 = facets[12*x+y,0:3]
                    c2 = facets[12*x+y,6:9]
                    c3 = facets[12*x+y,3:6]

                else:
                    InitialN = facetNormals[12*x+y,0:3]
                    c1 = facets[12*x+y,0:3]
                    c2 = facets[12*x+y,3:6]
                    c3 = facets[12*x+y,6:9]

                coords.append(np.array([c1,c2,c3]))
                InitialNormal.append(InitialN)


        InitialNormal=np.array(InitialNormal).reshape(-1,3)
        coords = np.array(coords).reshape(-1,9)
        facetNormals=InitialNormal
        facets = coords.reshape(-1,9)

        # Formation of rotation stacked matrix (3 x 3 x nFacets)
        v = np.cross(facetNormals,pn.reshape(-1,3))
        zeros = np.zeros(len(v),)
        ssc = np.array(([[zeros, -v[:,2], v[:,1]],[ v[:,2], zeros, -v[:,0]],\
            [ -v[:,1], v[:,0], zeros]]))
        identity = np.dstack([np.eye(3)]*len(v))
        mulNormalsPn = np.matmul(np.expand_dims(facetNormals, axis=2),\
            np.expand_dims(pn.reshape(-1,3), axis=1)).T
        R = identity + ssc + (np.matmul(ssc.T,ssc.T).T)*(1-mulNormalsPn)/\
        (np.dot(np.linalg.norm(v),np.linalg.norm(v)))

        # Clear not needed variables from memory
        del Check
        del c1
        del c2
        del c3
        del v
        del zeros
        del identity
        del mulNormalsPn
        del ssc        

        # Define 1st projected tangential
        tan1 = np.expand_dims(facetCenters-facets[:,0:3], axis=2)
        ptan1 = np.squeeze(np.matmul(np.transpose(R.T,(0, 2, 1)),tan1))/\
            np.array([np.linalg.norm(np.squeeze(np.matmul(np.transpose(R.T,\
                (0, 2, 1)),tan1)),axis=1),]*3).T

        # Define 2nd projected tangential
        ptan2 = np.cross(pn,ptan1)/np.array([np.linalg.norm(np.cross(pn,ptan1),\
            axis=1),]*3).T

        # Sub-tet Volume
        coord1 = facets[:,0:3]
        coord2 = facets[:,3:6]
        coord3 = facets[:,6:9]
        volCalc1 = np.expand_dims(np.transpose(coord1-p1).T, axis=1)
        volCalc2 = np.expand_dims(np.transpose(np.cross((coord2-p1),\
            (coord3-p1))).T, axis=2)
        volCalc3 = np.expand_dims(np.transpose(coord1-p2).T, axis=1)
        volCalc4 = np.expand_dims(np.transpose(np.cross((coord2-p2),\
            (coord3-p2))).T, axis=2)
        facetVol1 = np.squeeze(abs(np.matmul(volCalc1,volCalc2))/6)
        facetVol2 = np.squeeze(abs(np.matmul(volCalc3,volCalc4))/6)
        subtetVol = np.squeeze(abs(np.matmul(volCalc1,volCalc2))/6 \
            + abs(np.matmul(volCalc3,volCalc4))/6)
        #print(facetVol1)
        #print(facetVol2)

        # Clear not needed variables from memory
        del R
        del tan1
        del coord1
        del coord2
        del coord3
        del p1
        del p2
        del volCalc1
        del volCalc2
        del volCalc3
        del volCalc4

        # Projected Area
        areaCalc1 = np.squeeze(np.matmul(np.expand_dims(facetNormals, axis=1),\
            np.expand_dims(pn, axis=2)))
        areaCalc2 = np.linalg.norm(facetNormals, axis=1)\
            *np.linalg.norm(pn, axis=1)
        pArea = abs(areaCalc1/areaCalc2*facetAreas)

        # Initialize a data matrix for all facet data
        dataList = np.empty([len(allTets)*36,10])

        # Initialize a matrix for facet material information
        facetMaterial = np.empty([len(allTets)*12,])

        # Initialize a matrix for particle material information
        particleMaterial = np.empty([len(allTets)*12,2])


        # Extend material list for edge nodes
        if multiMaterial in ['on','On','Y','y','Yes','yes']:

            materialList = np.concatenate((2*np.ones([len(allNodes)-\
                len(materialList),]),materialList))

        elif cementStructure in ['on','On','Y','y','Yes','yes']:

        # materialList = np.concatenate((0*np.ones([len(allNodes)-\
        #     len(materialList),]),materialList))  
            materialList = np.concatenate((edgeMaterialList,materialList))

        else:

            materialList = np.concatenate((0*np.ones([len(allNodes)-\
                len(materialList),]),materialList))            

        for x in range(0,len(allTets)):

            for y in range(0,12):

                # Store: [n1 n2 cx cy cz fa nx ny nz vol]
                dataList[36*x+3*y,0] = allTets[x,tetn1[y]].astype(int)
                dataList[36*x+3*y,1] = allTets[x,tetn2[y]].astype(int)
                dataList[36*x+3*y,2:5] = facetCenters[12*x+y,:]
                dataList[36*x+3*y,5] = facetAreas[12*x+y]
                dataList[36*x+3*y,6:9] = facetNormals[12*x+y,:]
                dataList[36*x+3*y,9] = subtetVol[12*x+y]

                # Store: [pa px py pz qx qy qz sx sy sz]
                dataList[36*x+3*y+1,0] = pArea[12*x+y]
                dataList[36*x+3*y+1,1:4] = pn[12*x+y,:]
                dataList[36*x+3*y+1,4:7] = ptan1[12*x+y,:]
                dataList[36*x+3*y+1,7:10] = ptan2[12*x+y,:]

                # Store: [cx1 cy1 cz1 cx2 cy2 cz2 cx3 cy3 cz3 0]
            # if multiMaterial in ['on','On','Y','y','Yes','yes']
                # 1 = ITZ
                # 2 = Binder
                # 3 = Aggregate
                dataList[36*x+3*y+2,0:9] = facets[12*x+y,:]

                # Store particle materials
                particleMaterial[12*x+y,0] = materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)
                particleMaterial[12*x+y,1] = materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)

                # Material rule based on particle diameters volumes
                if materialRule == 9 or materialRule == 10:

                    if facetVol1.all() == facetVol2.all():
                        
                        dataList[36*x+3*y+2,9] = materialList[allTets[x,tetn1[y]]\
                            .astype(int)-1].astype(int)
                        facetMaterial[12*x+y] = materialList[allTets[x,tetn1[y]]\
                            .astype(int)-1].astype(int)

                    else:

                        if facetVol1.any() > facetVol2.any():

                            dataList[36*x+3*y+2,9] = materialList[allTets[x,tetn1[y]]\
                                .astype(int)-1].astype(int)
                            facetMaterial[12*x+y] = materialList[allTets[x,tetn1[y]]\
                                .astype(int)-1].astype(int)

                        else: 

                            dataList[36*x+3*y+2,9] = materialList[allTets[x,tetn2[y]]\
                                .astype(int)-1].astype(int)
                            facetMaterial[12*x+y] = materialList[allTets[x,tetn2[y]]\
                                .astype(int)-1].astype(int)

                elif materialRule > 0:

                    if materialRule == 1:
                        aggITZ = 3
                        aggBinder = 3
                        itzBinder = 1
                    elif materialRule == 2:
                        aggITZ = 3
                        aggBinder = 3
                        itzBinder = 2
                    elif materialRule == 3:
                        aggITZ = 3
                        aggBinder = 2
                        itzBinder = 1
                    elif materialRule == 4:
                        aggITZ = 3
                        aggBinder = 2
                        itzBinder = 2
                    elif materialRule == 5:
                        aggITZ = 1
                        aggBinder = 3
                        itzBinder = 1
                    elif materialRule == 6:
                        aggITZ = 1
                        aggBinder = 3
                        itzBinder = 2
                    elif materialRule == 7:
                        aggITZ = 1
                        aggBinder = 2
                        itzBinder = 1
                    elif materialRule == 8:
                        aggITZ = 1
                        aggBinder = 2
                        itzBinder = 2
                    
                    if materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int) \
                        == materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int):
                        
                        dataList[36*x+3*y+2,9] = materialList[allTets[x,tetn1[y]]\
                            .astype(int)-1].astype(int)
                        facetMaterial[12*x+y] = materialList[allTets[x,tetn1[y]]\
                            .astype(int)-1].astype(int)
                    
                    elif materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)==1 and materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)==2:

                        dataList[36*x+3*y+2,9] = itzBinder
                        facetMaterial[12*x+y] = itzBinder

                    elif materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)==2 and materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)==1:

                        dataList[36*x+3*y+2,9] = itzBinder
                        facetMaterial[12*x+y] = itzBinder

                    elif materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)==1 and materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)==3:

                        dataList[36*x+3*y+2,9] = aggITZ
                        facetMaterial[12*x+y] = aggITZ

                    elif materialList[allTets[x,tetn1[y]].astype(int)-1].astype(int)==3 and materialList[allTets[x,tetn2[y]].astype(int)-1].astype(int)==1:

                        dataList[36*x+3*y+2,9] = aggITZ
                        facetMaterial[12*x+y] = aggITZ

                    else:
                        
                        dataList[36*x+3*y+2,9] = aggBinder
                        facetMaterial[12*x+y] = aggBinder

                else:

                    dataList[36*x+3*y+2,9] = 0
                    facetMaterial[12*x+y] = 0                       

        return dataList,facetMaterial,subtetVol,facetVol1,facetVol2,particleMaterial


    def externalFacetFile(self,dataList,vertices,triangles,geoName):

        subtriangles = np.empty((len(triangles)*6,10))
        triangles = triangles.astype(int)

        for x in range(len(triangles)):

            # Face triangle
            #      n0
            #     / \
            #    /   \
            #  n1-----n2

            # [nx ny nz ex ey ez fx fy fz]
            # Node: n0
            subtriangles[6*x,:] = np.concatenate(([triangles[x,0]],vertices[triangles[x,0]-1],(vertices[triangles[x,0]-1]+vertices[triangles[x,1]-1])/2,((vertices[triangles[x,0]-1]+(vertices[triangles[x,1]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,1]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,2]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,1]-1])/2)/2)/3))
            subtriangles[6*x+1,:] = np.concatenate(([triangles[x,0]],vertices[triangles[x,0]-1],(vertices[triangles[x,0]-1]+vertices[triangles[x,2]-1])/2,((vertices[triangles[x,0]-1]+(vertices[triangles[x,1]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,1]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,2]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,1]-1])/2)/2)/3))

            # Node: n1
            subtriangles[6*x+2,:] = np.concatenate(([triangles[x,1]],vertices[triangles[x,1]-1],(vertices[triangles[x,0]-1]+vertices[triangles[x,1]-1])/2,((vertices[triangles[x,0]-1]+(vertices[triangles[x,1]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,1]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,2]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,1]-1])/2)/2)/3))
            subtriangles[6*x+3,:] = np.concatenate(([triangles[x,1]],vertices[triangles[x,1]-1],(vertices[triangles[x,1]-1]+vertices[triangles[x,2]-1])/2,((vertices[triangles[x,0]-1]+(vertices[triangles[x,1]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,1]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,2]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,1]-1])/2)/2)/3))

            # Node: n2
            subtriangles[6*x+4,:] = np.concatenate(([triangles[x,2]],vertices[triangles[x,2]-1],(vertices[triangles[x,0]-1]+vertices[triangles[x,2]-1])/2,((vertices[triangles[x,0]-1]+(vertices[triangles[x,1]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,1]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,2]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,1]-1])/2)/2)/3))
            subtriangles[6*x+5,:] = np.concatenate(([triangles[x,2]],vertices[triangles[x,2]-1],(vertices[triangles[x,1]-1]+vertices[triangles[x,2]-1])/2,((vertices[triangles[x,0]-1]+(vertices[triangles[x,1]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,1]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,2]-1])/2)/2+(vertices[triangles[x,2]-1]+(vertices[triangles[x,0]-1]+vertices[triangles[x,1]-1])/2)/2)/3))

        np.savetxt(Path('meshes/' + geoName + '/' + geoName \
            + '-data-faceFacets.dat'), subtriangles, fmt='%.10g', delimiter=' '\
            ,header='\
Facet Data Generated with LDPM Mesh Generation Tool\n\
Format: One line per sub triangle (surface facet)\n\
Note: First coordinate is corresponding tet node\n\
[n nx ny nz ex ey ez fx fy fz]')


    def materialSorting(self,facetMaterial,facetVol1,facetVol2,particleMaterial,subtetVol):
        
        # [#, Volume Difference, Selected Material, Material 1, Material 2, subtetVol]
        condensedData = np.vstack((np.arange(len(facetMaterial)),abs(facetVol1-facetVol2),facetMaterial,particleMaterial[:,0],particleMaterial[:,1],subtetVol)).transpose()

        sortedData = condensedData[condensedData[:,1].argsort()]

        return sortedData


    def materialRefinement(self,sortedData,itzVolFracSim,\
        binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct,i):
        
        # Same material case
        if sortedData[i,3] == sortedData[i,4]:
            #print('0')
            pass

        # Aggregate-ITZ Case
        elif (sortedData[i,3] == 1 and sortedData[i,4] == 3) or \
            (sortedData[i,3] == 3 and sortedData[i,4] == 1):

            if abs(itzVolFracSim-itzVolFracAct) > abs(aggVolFracSim-aggVolFracAct):
                #print('1')
                if itzVolFracSim-itzVolFracAct > 0:
                    sortedData[i,2] = 3
                else:
                    sortedData[i,2] = 1

            else:
                if aggVolFracSim-aggVolFracAct > 0:
                    sortedData[i,2] = 1
                else:
                    sortedData[i,2] = 3

        # Aggregate-Binder Case
        elif (sortedData[i,3] == 2 and sortedData[i,4] == 3) or \
            (sortedData[i,3] == 3 and sortedData[i,4] == 2):

            if abs(binderVolFracSim-binderVolFracAct) > abs(aggVolFracSim-aggVolFracAct):
                #print('2')
                if binderVolFracSim-binderVolFracAct > 0:
                    sortedData[i,2] = 3
                else:
                    sortedData[i,2] = 2

            else:
                if aggVolFracSim-aggVolFracAct > 0:
                    sortedData[i,2] = 2
                else:
                    sortedData[i,2] = 3

        # ITZ-Binder Case
        elif (sortedData[i,3] == 1 and sortedData[i,4] == 2) or \
            (sortedData[i,3] == 2 and sortedData[i,4] == 1):

            if abs(itzVolFracSim-itzVolFracAct) > abs(binderVolFracSim-binderVolFracAct):
                #print('3')
                if itzVolFracSim-itzVolFracAct > 0:
                    sortedData[i,2] = 2
                else:
                    sortedData[i,2] = 1

            else:
                if binderVolFracSim-binderVolFracAct > 0:
                    sortedData[i,2] = 1
                else:
                    sortedData[i,2] = 2
            
        return sortedData

    def cementMaterialRefinement(self,sortedData,PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,\
        CSH_LDVolFracSim,CSH_HDVolFracSim,PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,\
        CSH_LDVolFracAct,CSH_HDVolFracAct,i):
        
        # Same material case
        if sortedData[i,3] == sortedData[i,4]:
            #print('0')
            pass

        # Pores-Clinker Case
        elif (sortedData[i,3] == 1 and sortedData[i,4] == 2) or \
            (sortedData[i,3] == 2 and sortedData[i,4] == 1):

            if abs(PoresVolFracSim-PoresVolFracAct) > abs(ClinkerVolFracSim-ClinkerVolFracAct):
                #print('1')
                if PoresVolFracSim-PoresVolFracAct > 0:
                    sortedData[i,2] = 2
                else:
                    sortedData[i,2] = 1

            else:
                if ClinkerVolFracSim-ClinkerVolFracAct > 0:
                    sortedData[i,2] = 1
                else:
                    sortedData[i,2] = 2

        # Pores-CH Case
        elif (sortedData[i,3] == 1 and sortedData[i,4] == 3) or \
            (sortedData[i,3] == 3 and sortedData[i,4] == 1):

            if abs(PoresVolFracSim-PoresVolFracAct) > abs(CHVolFracSim-CHVolFracAct):
                #print('2')
                if PoresVolFracSim-PoresVolFracAct > 0:
                    sortedData[i,2] = 3
                else:
                    sortedData[i,2] = 1

            else:
                if CHVolFracSim-CHVolFracAct > 0:
                    sortedData[i,2] = 1
                else:
                    sortedData[i,2] = 3

        # Pores-CSH_LD Case
        elif (sortedData[i,3] == 1 and sortedData[i,4] == 4) or \
            (sortedData[i,3] == 4 and sortedData[i,4] == 1):

            if abs(PoresVolFracSim-PoresVolFracAct) > abs(CSH_LDVolFracSim-CSH_LDVolFracAct):
                #print('3')
                if PoresVolFracSim-PoresVolFracAct > 0:
                    sortedData[i,2] = 4
                else:
                    sortedData[i,2] = 1

            else:
                if CSH_LDVolFracSim-CSH_LDVolFracAct > 0:
                    sortedData[i,2] = 1
                else:
                    sortedData[i,2] = 4

        # Pores-CSH_HD Case
        elif (sortedData[i,3] == 1 and sortedData[i,4] == 5) or \
            (sortedData[i,3] == 5 and sortedData[i,4] == 1):

            if abs(PoresVolFracSim-PoresVolFracAct) > abs(CSH_HDVolFracSim-CSH_HDVolFracAct):
                #print('4')
                if PoresVolFracSim-PoresVolFracAct > 0:
                    sortedData[i,2] = 5
                else:
                    sortedData[i,2] = 1

            else:
                if CSH_HDVolFracSim-CSH_HDVolFracAct > 0:
                    sortedData[i,2] = 1
                else:
                    sortedData[i,2] = 5

        # Clinker-CH Case
        elif (sortedData[i,3] == 2 and sortedData[i,4] == 3) or \
            (sortedData[i,3] == 3 and sortedData[i,4] == 2):

            if abs(ClinkerVolFracSim-ClinkerVolFracAct) > abs(CHVolFracSim-CHVolFracAct):
                #print('5')
                if ClinkerVolFracSim-ClinkerVolFracAct > 0:
                    sortedData[i,2] = 3
                else:
                    sortedData[i,2] = 2

            else:
                if CHVolFracSim-CHVolFracAct > 0:
                    sortedData[i,2] = 2
                else:
                    sortedData[i,2] = 3

        # Clinker-CSH_LD Case
        elif (sortedData[i,3] == 2 and sortedData[i,4] == 4) or \
            (sortedData[i,3] == 4 and sortedData[i,4] == 2):

            if abs(ClinkerVolFracSim-ClinkerVolFracAct) > abs(CSH_LDVolFracSim-CSH_LDVolFracAct):
                #print('6')
                if ClinkerVolFracSim-ClinkerVolFracAct > 0:
                    sortedData[i,2] = 4
                else:
                    sortedData[i,2] = 2

            else:
                if CSH_LDVolFracSim-CSH_LDVolFracAct > 0:
                    sortedData[i,2] = 2
                else:
                    sortedData[i,2] = 4 

        # Clinker-CSH_HD Case
        elif (sortedData[i,3] == 2 and sortedData[i,4] == 5) or \
            (sortedData[i,3] == 5 and sortedData[i,4] == 2):

            if abs(ClinkerVolFracSim-ClinkerVolFracAct) > abs(CSH_HDVolFracSim-CSH_HDVolFracAct):
                #print('7')
                if ClinkerVolFracSim-ClinkerVolFracAct > 0:
                    sortedData[i,2] = 5
                else:
                    sortedData[i,2] = 2

            else:
                if CSH_HDVolFracSim-CSH_HDVolFracAct > 0:
                    sortedData[i,2] = 2
                else:
                    sortedData[i,2] = 5

        # CH-CSH_LD Case
        elif (sortedData[i,3] == 3 and sortedData[i,4] == 4) or \
            (sortedData[i,3] == 4 and sortedData[i,4] == 3):

            if abs(CHVolFracSim-CHVolFracAct) > abs(CSH_LDVolFracSim-CSH_LDVolFracAct):
                #print('8')
                if CHVolFracSim-CHVolFracAct > 0:
                    sortedData[i,2] = 4
                else:
                    sortedData[i,2] = 3

            else:
                if CSH_LDVolFracSim-CSH_LDVolFracAct > 0:
                    sortedData[i,2] = 3
                else:
                    sortedData[i,2] = 4   

        # CH-CSH_HD Case
        elif (sortedData[i,3] == 3 and sortedData[i,4] == 5) or \
            (sortedData[i,3] == 5 and sortedData[i,4] == 3):

            if abs(CHVolFracSim-CHVolFracAct) > abs(CSH_HDVolFracSim-CSH_HDVolFracAct):
                #print('9')
                if CHVolFracSim-CHVolFracAct > 0:
                    sortedData[i,2] = 5
                else:
                    sortedData[i,2] = 3

            else:
                if CSH_HDVolFracSim-CSH_HDVolFracAct > 0:
                    sortedData[i,2] = 3
                else:
                    sortedData[i,2] = 5

        # CSH_LD-CSH_HD Case
        elif (sortedData[i,3] == 4 and sortedData[i,4] == 5) or \
            (sortedData[i,3] == 5 and sortedData[i,4] == 4):

            if abs(CSH_LDVolFracSim-CSH_LDVolFracAct) > abs(CSH_HDVolFracSim-CSH_HDVolFracAct):
                #print('10')
                if CSH_LDVolFracSim-CSH_LDVolFracAct > 0:
                    sortedData[i,2] = 5
                else:
                    sortedData[i,2] = 4

            else:
                if CSH_HDVolFracSim-CSH_HDVolFracAct > 0:
                    sortedData[i,2] = 4
                else:
                    sortedData[i,2] = 5

        return sortedData

    def reformDataList(self,allTets,dataList,sortedData):

        reSortedData = sortedData[sortedData[:,0].argsort()]

        facetMaterial = reSortedData[:,2]

        for x in range(0,len(allTets)):

            for y in range(0,12):

                dataList[36*x+3*y+2,9] = facetMaterial[12*x+y]
                
        return dataList,facetMaterial


    def facetFile(self,geoName,dataList,allTets,dataType):
        

        if dataType in ['t','T','text','Text','raw','Raw']:

            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-facet.dat'), dataList, fmt='%.10g', delimiter=' '\
                ,header='\
Facet Data Generated with LDPM Mesh Generation Tool\n\
Format: Tet label followed by three lines per facet\n\
Tet [N]\n\
[n1 n2 cx cy cz fa nx ny nz vol]\n\
[pa px py pz qx qy qz sx sy sz]\n\
[fx1 fy1 fz1 fx2 fy2 fz2 fx3 fy3 fz3 mF]')

            f = open(Path('meshes/' + geoName + '/' + geoName \
                + '-data-facet.dat'), "r")
            contents = f.readlines()
            f.close()

            contents.insert(6, '\n')
            contents.insert(7, 'Number of LDPM tets: '+str(len(allTets)) + '\n')
            for x in range(0,len(allTets)):
                contents.insert(37*x+8, 'Tet index: ' + str(x+1) + '\n')

            f = open(Path('meshes/' + geoName + '/' + geoName \
                + '-data-facet.dat'), "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()

        else:
            np.save(Path('meshes/' + geoName + '/' + geoName \
                + '-data-facet.dat'), dataList, allow_pickle=False)


    def facetRandomFieldFile(self,geoName,dataList,allTets,dataType):

        RandFile = []
        for x in range(0,len(allTets)):
            for y in range(0,12):
                RandFile.append(np.array([dataList[12*x+y,2],dataList[12*x+y,3],dataList[12*x+y,4],int(x+1),int(y+1)]))

        RandFile = np.array(RandFile).reshape(-1,5)
        
        if dataType in ['t','T','text','Text','raw','Raw']:

            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-facetRandomfieldFile.dat'), RandFile, fmt='%.10g', delimiter=' ')

        else:
            np.save(Path('meshes/' + geoName + '/' + geoName \
                + '-data-facetRandomfieldFile.dat'), RandFile, allow_pickle=False)


    def facetVTKFile(self,geoName,tetFacets,dataList,facetMaterial,\
        multiMaterial,cementStructure):

        facetsPoints = tetFacets.reshape(-1,3)
        cells = (np.arange(0,round(len(facetsPoints))).\
            reshape(-1,3)).astype(int)
        cell_types = np.array([5,]*round(len(facetsPoints)/3))

        with open(Path('meshes/' + geoName + '/' + geoName + \
            '-para-facet.000.vtk'),"w") as f:                                                                          
            f.write('# vtk DataFile Version 2.0\n')
            f.write('Unstructured Grid\n')            
            f.write('ASCII\n')    
            f.write('DATASET UNSTRUCTURED_GRID\n')        
            f.write('POINTS ' + str(len(facetsPoints)) + ' double \n')  
            f.write("\n".join(" ".join(map(str, x)) for x in facetsPoints))
            f.write('\n\n')  
            f.write('CELLS ' + str(round(len(facetsPoints)/3)) + ' ' \
                + str(round(len(facetsPoints)/3*4)) +'\n3 ')
            f.write("\n3 ".join(" ".join(map(str, x)) for x in cells))
            f.write('\n\n')  
            f.write('CELL_TYPES ' + str(round(len(facetsPoints)/3)) +'\n')
            for x in cell_types:
                f.write("%s\n" % x)
            if multiMaterial in ['on','On','Y','y','Yes','yes']:  
                f.write('\nCELL_DATA ' + str(len(facetMaterial)) + '\n')
                f.write('FIELD FieldData 1\n')
                f.write('material 1 ' + str(len(facetMaterial)) + ' float\n')
                for x in facetMaterial:
                    f.write("%s\n" % x)
            if cementStructure in ['on','On','Y','y','Yes','yes']:  
                f.write('\nCELL_DATA ' + str(len(facetMaterial)) + '\n')
                f.write('FIELD FieldData 1\n')
                f.write('material 1 ' + str(len(facetMaterial)) + ' float\n')
                for x in facetMaterial:
                    f.write("%s\n" % x) 

    def facetfiberInteractionData(self,p1Fiber,p2Fiber,dFiber,lFiber,orienFibers,\
        geoName,allTets,allNodes,tetFacets,dataList,tetn1,tetn2,facetNormals,facetCenters):
        # Number of total fiber
        NumberofFibers = len(p1Fiber[:,1])
        print('NumberofFibers',NumberofFibers)
        # Initialize a data array and matrix for interseted fiber
        FiberdataList = []
        ProjectedFacet = []
        FibertetList = []
        IntersectedFiber = []
        No = 0.0
        TotalTet = 0.0
        TotalFiber = 0.0


        for z in range(0,NumberofFibers):
            # Obtain extents for floating bin for fiber
            binMin = np.amin(np.vstack((p1Fiber[z,:],p2Fiber[z,:])), axis=0)-1.0*lFiber[z]-dFiber
            binMax = np.amax(np.vstack((p1Fiber[z,:],p2Fiber[z,:])), axis=0)+1.0*lFiber[z]+dFiber

            for x in range(0,len(allTets)):
                # Store tet vertices that fall inside the bin

                coord1 = allNodes[(allTets[x,0]-1).astype(int)]
                coord2 = allNodes[(allTets[x,1]-1).astype(int)]
                coord3 = allNodes[(allTets[x,2]-1).astype(int)]
                coord4 = allNodes[(allTets[x,3]-1).astype(int)]

                coord1 = np.all([(coord1[0] > binMin[0]) , (coord1[0] < binMax[0]),\
                    (coord1[1] > binMin[1]) , (coord1[1] < binMax[1]) ,\
                    (coord1[2] > binMin[2]) , (coord1[2] < binMax[2])],axis=0)      
                coord2 = np.all([(coord2[0] > binMin[0]) , (coord2[0] < binMax[0]),\
                    (coord2[1] > binMin[1]) , (coord2[1] < binMax[1]) ,\
                    (coord2[2] > binMin[2]) , (coord2[2] < binMax[2])],axis=0)          
                coord3 = np.all([(coord3[0] > binMin[0]) , (coord3[0] < binMax[0]),\
                    (coord3[1] > binMin[1]) , (coord3[1] < binMax[1]) ,\
                    (coord3[2] > binMin[2]) , (coord3[2] < binMax[2])],axis=0)  
                coord4 = np.all([(coord4[0] > binMin[0]) , (coord4[0] < binMax[0]),\
                    (coord4[1] > binMin[1]) , (coord4[1] < binMax[1]) ,\
                    (coord4[2] > binMin[2]) , (coord4[2] < binMax[2])],axis=0)  

                binTets = np.any([coord1,coord2,coord3,coord4],axis=0) 
                
                if binTets == True:
                    FibertetList.append(np.array([int(x),int(z)]))  

        FibertetList=np.array(FibertetList).reshape(-1,2)
        FibertetList=FibertetList[FibertetList[:,0].argsort()]                          
        FiberBin=np.unique(FibertetList[:,1])
        print('FiberBin',len(FiberBin))

        #######################################################
        # For fiberfacet using Matthew's rotational matrix

        # Store facets as set of three points
        facets = tetFacets.reshape(-1,9)

        # Store particles accociated with facets
        p1 = allNodes[(allTets[:,tetn1]-1).astype(int),:].reshape(-1,3)
        p2 = allNodes[(allTets[:,tetn2]-1).astype(int),:].reshape(-1,3)


        # Projected facet normal
        pn = (p2-p1)/np.array([np.linalg.norm(p2-p1,axis=1),]*3).T
        pn = pn.reshape(-1,3)
        #print('n',pn)

        InitialNormal=[]
        coords=[]

        for x in range(0,len(allTets)):

            for y in range(0,12):


                Check = pn[12*x+y,0]*facetNormals[12*x+y,0]+\
                        pn[12*x+y,1]*facetNormals[12*x+y,1]+\
                        pn[12*x+y,2]*facetNormals[12*x+y,2]

                if Check < 0.0:
                    InitialN = -1.0*facetNormals[12*x+y,0:3]
                    c1 = facets[12*x+y,0:3]
                    c2 = facets[12*x+y,6:9]
                    c3 = facets[12*x+y,3:6]

                else:
                    InitialN = facetNormals[12*x+y,0:3]
                    c1 = facets[12*x+y,0:3]
                    c2 = facets[12*x+y,3:6]
                    c3 = facets[12*x+y,6:9]

                coords.append(np.array([c1,c2,c3]))
                InitialNormal.append(InitialN)

        InitialNormal = np.array(InitialNormal).reshape(-1,3)
        coords = np.array(coords).reshape(-1,9)
        facetNormals = InitialNormal
        facets = coords.reshape(-1,9)

        # Formation of rotation stacked matrix (3 x 3 x nFacets)
        v = np.cross(facetNormals,pn.reshape(-1,3))
        zeros = np.zeros(len(v),)
        ssc = np.array(([[zeros, -v[:,2], v[:,1]],[ v[:,2], zeros, -v[:,0]],\
            [ -v[:,1], v[:,0], zeros]]))
        identity = np.dstack([np.eye(3)]*len(v))
        mulNormalsPn = np.matmul(np.expand_dims(facetNormals, axis=2),\
            np.expand_dims(pn.reshape(-1,3), axis=1)).T
        R = identity + ssc + (np.matmul(ssc.T,ssc.T).T)*(1-mulNormalsPn)/\
        (np.dot(np.linalg.norm(v),np.linalg.norm(v)))

        # Make vectors from center to corners of the facets
        vectorOA = np.expand_dims(facetCenters[:,0:3]-facets[:,0:3], axis=2)
        vectorOB = np.expand_dims(facetCenters[:,0:3]-facets[:,3:6], axis=2)
        vectorOC = np.expand_dims(facetCenters[:,0:3]-facets[:,6:9], axis=2)

        pvectorOA =np.squeeze(np.matmul(np.transpose(R.T,(0, 2, 1)),vectorOA))
        pvectorOB = np.squeeze(np.matmul(np.transpose(R.T,(0, 2, 1)),vectorOB))
        pvectorOC = np.squeeze(np.matmul(np.transpose(R.T,(0, 2, 1)),vectorOC))

        RcoordP1= facetCenters-pvectorOA
        RcoordP2= facetCenters-pvectorOB
        RcoordP3= facetCenters-pvectorOC

        # Compute the plane supporting the triangle (coordP1, coordP2, coordP3)  normal: nprojectedfacet offset: d
        vector12 = RcoordP2-RcoordP1
        vector23 = RcoordP3-RcoordP2
        vector13 = RcoordP3-RcoordP1
        vector31 = RcoordP1-RcoordP3

        Normal = np.cross(vector12,vector13)/np.array([np.linalg.norm(np.cross(vector12,vector13),axis=1),]*3).T

        RcoordP1=np.array(RcoordP1).reshape(-1,3)
        RcoordP2=np.array(RcoordP2).reshape(-1,3)
        RcoordP3=np.array(RcoordP3).reshape(-1,3)
        vector12=np.array(vector12).reshape(-1,3)
        vector31=np.array(vector31).reshape(-1,3)
        vector23=np.array(vector23).reshape(-1,3)
        vector13=np.array(vector13).reshape(-1,3)
        Normal=np.array(Normal).reshape(-1,3)
        projectedFacet = np.concatenate((RcoordP1,RcoordP2,RcoordP3))

        for i in range(0,len(FibertetList)):

            x=int(FibertetList[i,0])
            z=int(FibertetList[i,1])

            TotalTet = TotalTet + 1    

            # Check the intersection of this inside fiber with 12 facsets of tet
            for y in range(0,12):

                tan1 = dataList[36*x+3*y+1,4:7]
                tan2 = dataList[36*x+3*y+1,7:10]

                p12Fiber = np.array([p2Fiber[z,0]-p1Fiber[z,0],\
                    p2Fiber[z,1]-p1Fiber[z,1],\
                    p2Fiber[z,2]-p1Fiber[z,2]])


                offsetd = Normal[12*x+y,0]*RcoordP1[12*x+y,0]+\
                    Normal[12*x+y,1]*RcoordP1[12*x+y,1]+\
                    Normal[12*x+y,2]*RcoordP1[12*x+y,2]

                offsetd = -1*offsetd

                facetnormalDotp12Fiber = Normal[12*x+y,0]*p12Fiber[0]+\
                    Normal[12*x+y,1]*p12Fiber[1]+\
                    Normal[12*x+y,2]*p12Fiber[2]
                

                facetnormalDotp1Fiber = Normal[12*x+y,0]*p1Fiber[z,0]+\
                    Normal[12*x+y,1]*p1Fiber[z,1]+\
                    Normal[12*x+y,2]*p1Fiber[z,2]

                # Ignore line parallel to (or lying in) the plane
                if abs(facetnormalDotp12Fiber) > 0.0:
                    t = - (offsetd + facetnormalDotp1Fiber)/facetnormalDotp12Fiber

                # Check if the intersection point is between p1Fiber and p2Fiber
                    if t >= 0.0 and t <= 1.0:

                        # From here can find why we have loop over tets, facets and fibers
                        P = np.array([p1Fiber[z,0] + t*p12Fiber[0],\
                            p1Fiber[z,1] + t*p12Fiber[1],\
                            p1Fiber[z,2] + t*p12Fiber[2]])

                        allPoints =[]
                        CheckEdge1 = []
                        CheckEdge2 = []
                        CheckEdge3 = []

                        # N12 = np.cross(Normal[12*x+y,0:3],vector12[12*x+y,0:3])
                        # PP1 = np.array([P[0]-RcoordP1[12*x+y,0],P[1]-RcoordP1[12*x+y,1],P[2]-RcoordP1[12*x+y,2]])
                        # CheckEdge1 = np.dot(PP1.T,N12)/np.linalg.norm(N12)

                        PP1 = np.array([P[0]-RcoordP1[12*x+y,0],P[1]-RcoordP1[12*x+y,1],P[2]-RcoordP1[12*x+y,2]])
                        N12 = np.cross(vector12[12*x+y,0:3],PP1)
                        CheckEdge1 = np.dot(Normal[12*x+y,0:3].T,N12)

                        if (CheckEdge1>=0).any():

                            # N23 = np.cross(Normal[12*x+y,0:3],vector23[12*x+y,0:3])
                            # PP2 = np.array([P[0]-RcoordP2[12*x+y,0],P[1]-RcoordP2[12*x+y,1],P[2]-RcoordP2[12*x+y,2]])
                            # CheckEdge2 = np.dot(PP2.T,N23)/np.linalg.norm(N23)

                            PP2 = np.array([P[0]-RcoordP2[12*x+y,0],P[1]-RcoordP2[12*x+y,1],P[2]-RcoordP2[12*x+y,2]])
                            N23 = np.cross(vector23[12*x+y,0:3],PP2)
                            CheckEdge2 = np.dot(Normal[12*x+y,0:3].T,N23)

                            if (CheckEdge2>=0).any():

                                # N31 = np.cross(Normal[12*x+y,0:3],vector31[12*x+y,0:3])
                                # PP3 = np.array([P[0]-RcoordP3[12*x+y,0],P[1]-RcoordP3[12*x+y,1],P[2]-RcoordP3[12*x+y,2]])
                                # CheckEdge3 = np.dot(PP3.T,N31)/np.linalg.norm(N31)

                                PP3 = np.array([P[0]-RcoordP3[12*x+y,0],P[1]-RcoordP3[12*x+y,1],P[2]-RcoordP3[12*x+y,2]])
                                N31 = np.cross(vector31[12*x+y,0:3],PP3)
                                CheckEdge3 = np.dot(Normal[12*x+y,0:3].T,N31)

                                if (CheckEdge3>=0).any():

                                   # Determination of Short and Long lenght of fiber intersected facet
                                    IntersectedFiber.append(np.array([int(z)]))
                                    No = No + 1
                                    distancep1Fiber = math.sqrt(((P[0]-p1Fiber[z,0])**2)+\
                                        ((P[1]-p1Fiber[z,1])**2)+\
                                        ((P[2]-p1Fiber[z,2])**2))
                                    distancep2Fiber = math.sqrt(((P[0]-p2Fiber[z,0])**2)+\
                                        ((P[1]-p2Fiber[z,1])**2)+\
                                        ((P[2]-p2Fiber[z,2])**2))

                                    tetindex = x+1
                                    facetindex = y+1
                                    FiberShortLenght = distancep2Fiber
                                    FiberLongLenght = distancep1Fiber

                                    if  distancep1Fiber < distancep2Fiber:
                                        FiberShortLenght = distancep1Fiber
                                        FiberLongLenght = distancep2Fiber

                                    InterPerFacet = 1.0

                                    OneIntersectedFiber = np.array([tetindex, facetindex, orienFibers[z,0], orienFibers[z,1], orienFibers[z,2],\
                                        FiberShortLenght, FiberLongLenght, dFiber, InterPerFacet])

                                    FiberdataList.append(OneIntersectedFiber)
                TotalIntersections = No

        FiberdataList = np.array(FiberdataList)
        FiberdataList = FiberdataList.reshape(-1,9)
        FiberdataList=FiberdataList[FiberdataList[:,0].argsort()]
        FiberdataList=FiberdataList[np.lexsort(FiberdataList[:,::-1].T)]
        IntersectedFiber = np.unique(np.array(IntersectedFiber))
        TotalFiber = len(IntersectedFiber)
        print('TotalFiber',TotalFiber)
        print(IntersectedFiber)


        for k in range(1,len(FiberdataList)):
            if ((FiberdataList[k,0]==FiberdataList[k-1,0]) and (FiberdataList[k,1]==FiberdataList[k-1,1])):
                FiberdataList[k,8] = FiberdataList[k-1,8] + 1

        MaxInterPerFacet = np.max(FiberdataList[:,8])

        return FiberdataList,TotalIntersections,MaxInterPerFacet,TotalTet,TotalFiber,IntersectedFiber,projectedFacet


    def facetfiberInteractionFile(self,geoName,FiberdataList,TotalIntersections,MaxInterPerFacet,TotalTet,TotalFiber,dataType):
    # Generate file for fiber-facet interaction data
        if dataType in ['t','T','text','Text','raw','Raw']:
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-fiberfacet.dat'), FiberdataList, fmt='%.10g', delimiter=' '\
                ,header='\
FacetFiber Data Generated with LDPM Mesh Generation Tool\n\
Format: Intersection number followed by one line per Interaction\n\
Intersection [N]\n\
[Te Fa fN fM fL S L df I]')

            f = open(Path('meshes/' + geoName + '/' + geoName \
                + '-data-fiberfacet.dat'), "r")
            contents = f.readlines()
            f.close()

            contents.insert(4, '\n')
            contents.insert(5, 'Number of Total Intersection: ' + str(int(TotalIntersections)) + '\n')
            contents.insert(6, 'Number of Maximum Intersection: ' + str(int(MaxInterPerFacet)) + '\n')
            for x in range(0,int(TotalIntersections)):
                contents.insert(2*x+7, 'Intersection: ' + str(x+1) + '\n')

            f = open(Path('meshes/' + geoName + '/' + geoName \
                + '-data-fiberfacet.dat'), "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()

        else:
            np.save(Path('meshes/' + geoName + '/' + geoName \
                + '-data-fiberfacet.dat'), dataList, allow_pickle=False)
            
    def sieveFile(self,geoName,multiMaterial,aggDiameterList,maxAggD,minAggD,aggGrainsDiameterList,\
        itzDiameterList,binderDiameterList,cementStructure,PoresDiameterList,ClinkerDiameterList,CHDiameterList,\
        CSH_LDDiameterList,CSH_HDDiameterList):
        
        if multiMaterial in ['on','On','Y','y','Yes','yes']: 

            # Generate sieve curve data file for aggregate
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurveAgg.dat'), aggGrainsDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):')           

            # Generate sieve curve data file for ITZ
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurveITZ.dat'), itzDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):')       


            # Generate sieve curve data file for binder
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurveBinder.dat'), binderDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):')   


            # Calculations for sieve curve plotting aggregate
            totalVol = sum(4/3*math.pi*(aggGrainsDiameterList/2)**3)

            passingAgg = np.zeros(len(aggGrainsDiameterList))

            for x in range(len(aggGrainsDiameterList)):
                passingAgg[x] = sum(4/3*math.pi*(aggGrainsDiameterList[0:(len(aggGrainsDiameterList)-x)]/2)**3)/totalVol*100

            # Calculations for sieve curve plotting ITZ
            totalVol = sum(4/3*math.pi*(itzDiameterList/2)**3)

            passingITZ = np.zeros(len(itzDiameterList))

            for x in range(len(itzDiameterList)):
                passingITZ[x] = sum(4/3*math.pi*(itzDiameterList[0:(len(itzDiameterList)-x)]/2)**3)/totalVol*100

            # Calculations for sieve curve plotting binder
            totalVol = sum(4/3*math.pi*(binderDiameterList/2)**3)

            passingBinder = np.zeros(len(binderDiameterList))

            for x in range(len(binderDiameterList)):
                passingBinder[x] = sum(4/3*math.pi*(binderDiameterList[0:(len(binderDiameterList)-x)]/2)**3)/totalVol*100

            # Generate plot of sieve curves
            plt.close('all')
            plt.plot(aggGrainsDiameterList, passingAgg, label= "Aggregate", color='k',linestyle='-') 
            plt.plot(itzDiameterList, passingITZ, label= "ITZ", color='0.8',linestyle='--') 
            plt.plot(binderDiameterList, passingBinder, label= "Binder", color='0.6',linestyle='-.') 


            font = FontProperties()
            font.set_family('serif')
            font.set_name('Times New Roman')
            plt.rcParams["mathtext.fontset"] = "dejavuserif"
            plt.rcParams["font.family"] = "Times New Roman"

            plt.title("Particle Sieve Curve", fontproperties=font)
            plt.xlabel('Particle Diameter, $d$ (mm)', fontproperties=font) 
            plt.ylabel('Percent Passing, $P$ (%)', fontproperties=font)

            plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
            plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))

            label_format = '{:,.1f}'
            ticks_loc = plt.gca().get_yticks().tolist()
            plt.gca().yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            plt.gca().set_yticklabels([label_format.format(x) for x in ticks_loc], fontproperties=font)

            ticks_loc = plt.gca().get_xticks().tolist()
            plt.gca().xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            plt.gca().set_xticklabels([label_format.format(x) for x in ticks_loc], fontproperties=font)


            plt.tick_params(axis='both', which='minor', labelsize=10)
            plt.tick_params(axis='both', which='major', labelsize=10)

            plt.grid(color='0.75', linestyle='--', linewidth=0.5)
            plt.legend() 

            plt.savefig(Path('meshes/' + geoName + '/' + geoName \
                + '-sieveCurve.png'), dpi=300)

        elif cementStructure in ['on','On','Y','y','Yes','yes']: 

            # Generate sieve curve data file for Capillary Pores
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurvePores.dat'), PoresDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):')           

            # Generate sieve curve data file for Clinker
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurveClinker.dat'), ClinkerDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):')       


            # Generate sieve curve data file for CH
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurveCH.dat'), CHDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):')   


            # Generate sieve curve data file for CSH_LD
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurveCSH_LD.dat'), CSH_LDDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):')   


            # Generate sieve curve data file for CSH_HD
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurveCSH_HD.dat'), CSH_HDDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):') 


            # Calculations for sieve curve plotting Capillary Pores
            totalVol = sum(4/3*math.pi*(PoresDiameterList/2)**3)

            passingPores = np.zeros(len(PoresDiameterList))

            for x in range(len(PoresDiameterList)):
                passingPores[x] = sum(4/3*math.pi*(PoresDiameterList[0:(len(PoresDiameterList)-x)]/2)**3)/totalVol*100

            # Calculations for sieve curve plotting Clinker
            totalVol = sum(4/3*math.pi*(ClinkerDiameterList/2)**3)

            passingClinker = np.zeros(len(ClinkerDiameterList))

            for x in range(len(ClinkerDiameterList)):
                passingClinker[x] = sum(4/3*math.pi*(ClinkerDiameterList[0:(len(ClinkerDiameterList)-x)]/2)**3)/totalVol*100

            # Calculations for sieve curve plotting CH
            totalVol = sum(4/3*math.pi*(CHDiameterList/2)**3)

            passingCH = np.zeros(len(CHDiameterList))

            for x in range(len(CHDiameterList)):
                passingCH[x] = sum(4/3*math.pi*(CHDiameterList[0:(len(CHDiameterList)-x)]/2)**3)/totalVol*100

            # Calculations for sieve curve plotting CSH_LD
            totalVol = sum(4/3*math.pi*(CSH_LDDiameterList/2)**3)

            passingCSH_LD = np.zeros(len(CSH_LDDiameterList))

            for x in range(len(CSH_LDDiameterList)):
                passingCSH_LD[x] = sum(4/3*math.pi*(CSH_LDDiameterList[0:(len(CSH_LDDiameterList)-x)]/2)**3)/totalVol*100

            # Calculations for sieve curve plotting CSH_HD
            totalVol = sum(4/3*math.pi*(CSH_HDDiameterList/2)**3)

            passingCSH_HD = np.zeros(len(CSH_HDDiameterList))

            for x in range(len(CSH_HDDiameterList)):
                passingCSH_HD[x] = sum(4/3*math.pi*(CSH_HDDiameterList[0:(len(CSH_HDDiameterList)-x)]/2)**3)/totalVol*100


            # Generate plot of sieve curves
            plt.close('all')
            plt.plot(PoresDiameterList, passingPores, label= "Pores", color='k',linestyle='--') 
            plt.plot(ClinkerDiameterList, passingClinker, label= "Clinker", color='0.8',linestyle='--') 
            plt.plot(CHDiameterList, passingCH, label= "CH", color='0.6',linestyle='--') 
            plt.plot(CSH_LDDiameterList, passingCSH_LD, label= "CSH_LD", color='0.4',linestyle='--') 
            plt.plot(CSH_HDDiameterList, passingCSH_HD, label= "CSH_HD", color='0.2',linestyle='--') 

            font = FontProperties()
            font.set_family('serif')
            font.set_name('Times New Roman')
            plt.rcParams["mathtext.fontset"] = "dejavuserif"
            plt.rcParams["font.family"] = "Times New Roman"

            plt.title("Particle Sieve Curve", fontproperties=font)
            plt.xlabel('Particle Diameter, $d$ (mm)', fontproperties=font) 
            plt.ylabel('Percent Passing, $P$ (%)', fontproperties=font)

            plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
            plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
            plt.gca().set_xticklabels(plt.gca().get_xticks(), fontproperties=font)
            plt.gca().set_yticklabels(plt.gca().get_yticks(), fontproperties=font)

            plt.tick_params(axis='both', which='minor', labelsize=10)
            plt.tick_params(axis='both', which='major', labelsize=10)

            plt.grid(color='0.75', linestyle='--', linewidth=0.5)
            plt.legend() 

            plt.savefig(Path('meshes/' + geoName + '/' + geoName \
                + '-sieveCurve.png'), dpi=300)
        else:

            # Generate sieve curve data file
            np.savetxt(Path('meshes/' + geoName + '/' + geoName \
                + '-data-sieveCurve.dat'), aggDiameterList, fmt='%.10g', delimiter=' '\
                ,header='\
    Data Generated with LDPM Mesh Generation Tool\n\
    Diameter List (mm):')   

            # Calculations for sieve curve plotting
            totalVol = sum(4/3*math.pi*(aggDiameterList/2)**3)

            passing = np.zeros(len(aggDiameterList))

            for x in range(len(aggDiameterList)):
                passing[x] = sum(4/3*math.pi*(aggDiameterList[0:(len(aggDiameterList)-x)]/2)**3)/totalVol*100
            

            # Generate plot of sieve curve
            plt.close('all')
            plt.plot(aggDiameterList, passing, label= "Simulated Data", color='k') 

            font = FontProperties()
            font.set_family('serif')
            font.set_name('Times New Roman')
            plt.rcParams["mathtext.fontset"] = "dejavuserif"
            plt.rcParams["font.family"] = "Times New Roman"

            plt.title("Particle Sieve Curve", fontproperties=font)

            plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
            plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))

            label_format = '{:,.1f}'
            ticks_loc = plt.gca().get_yticks().tolist()
            plt.gca().yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            plt.gca().set_yticklabels([label_format.format(x) for x in ticks_loc], fontproperties=font)

            ticks_loc = plt.gca().get_xticks().tolist()
            plt.gca().xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            plt.gca().set_xticklabels([label_format.format(x) for x in ticks_loc], fontproperties=font)

            plt.xlabel('Particle Diameter, $d$ (mm)', fontproperties=font) 
            plt.ylabel('Percent Passing, $P$ (%)', fontproperties=font)

            plt.tick_params(axis='both', which='minor', labelsize=10)
            plt.tick_params(axis='both', which='major', labelsize=10)

            plt.grid(color='0.75', linestyle='--', linewidth=0.5)
            plt.legend() 

            plt.savefig(Path('meshes/' + geoName + '/' + geoName \
                + '-sieveCurve.png'), dpi=300)

    def materialVolCheck(self,geoName,subtetVol,facetMaterial,agg,itz,binder):

        # Facet volumes
        itzVol = sum(subtetVol*(facetMaterial==1))

        binderVol = sum(subtetVol*(facetMaterial==2))

        aggVol = sum(subtetVol*(facetMaterial==3))

        # Facet volume fractions
        itzVolFracSim = itzVol/(sum((itzVol,binderVol,aggVol)))

        binderVolFracSim = binderVol/(sum((itzVol,binderVol,aggVol)))

        aggVolFracSim = aggVol/(sum((itzVol,binderVol,aggVol)))

        # Voxel volume fractions
        itzVolFracAct = len(itz)/(len(agg)+len(itz)+len(binder))

        binderVolFracAct = len(binder)/(len(agg)+len(itz)+len(binder))

        aggVolFracAct = len(agg)/(len(agg)+len(itz)+len(binder))

        return itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct

    def cementMaterialVolCheck(self,geoName,subtetVol,facetMaterial,PoresVoxelData,ClinkerVoxelData,\
        CHVoxelData,CSH_LDVoxelData,CSH_HDVoxelData,AllHydrationProduct):

        # Facet volumes
        PoresVol = sum(subtetVol*(facetMaterial==1))

        ClinkerVol = sum(subtetVol*(facetMaterial==2))

        CHVol = sum(subtetVol*(facetMaterial==3))

        CSH_LDVol = sum(subtetVol*(facetMaterial==4))

        CSH_HDVol = sum(subtetVol*(facetMaterial==5))

        # Facet volume fractions
        PoresVolFracSim = PoresVol/(sum((PoresVol,ClinkerVol,CHVol,CSH_LDVol,CSH_HDVol)))

        ClinkerVolFracSim = ClinkerVol/(sum((PoresVol,ClinkerVol,CHVol,CSH_LDVol,CSH_HDVol)))

        CHVolFracSim = CHVol/(sum((PoresVol,ClinkerVol,CHVol,CSH_LDVol,CSH_HDVol)))

        CSH_LDVolFracSim = CSH_LDVol/(sum((PoresVol,ClinkerVol,CHVol,CSH_LDVol,CSH_HDVol)))

        CSH_HDVolFracSim = CSH_HDVol/(sum((PoresVol,ClinkerVol,CHVol,CSH_LDVol,CSH_HDVol)))

        # Voxel volume fractions
        PoresVolFracAct = len(PoresVoxelData)/len(AllHydrationProduct)

        ClinkerVolFracAct = len(ClinkerVoxelData)/len(AllHydrationProduct)

        CHVolFracAct = len(CHVoxelData)/len(AllHydrationProduct)

        CSH_LDVolFracAct = len(CSH_LDVoxelData)/len(AllHydrationProduct)

        CSH_HDVolFracAct = len(CSH_HDVoxelData)/len(AllHydrationProduct)

        return PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
               PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct

    def materialVolPlot(self,geoName,itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,aggVolFracAct):

        # Plot volume fraction comparison
        plt.close('all')
        fig, ax = plt.subplots()

        labels = ['ITZ', 'Binder', 'Aggregate']

        simulated = [itzVolFracSim,binderVolFracSim,aggVolFracSim]

        voxelated = [itzVolFracAct,binderVolFracAct,aggVolFracAct]

        # Set up font
        font = FontProperties()
        font.set_family('serif')
        font.set_name('Times New Roman')
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "dejavuserif"

        # Set width of bar
        barWidth = 0.35
          
        # Set position of bar on X axis
        r1 = np.arange(len(simulated))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
         
        # Make the plot
        plt.bar(r1, simulated, color='0.6', width=barWidth, edgecolor='white', label='Simulated')
        plt.bar(r2, voxelated, color='0.8', width=barWidth, edgecolor='white', label='Voxelated')
        for i, v in enumerate(simulated):
            plt.text(i-0.07, v+0.01, "{:.2f}".format(v), color='k', fontproperties=font)
        for i, v in enumerate(voxelated):
            plt.text(i+0.28, v+0.01, "{:.2f}".format(v), color='k', fontproperties=font)
         
        # Add xticks on the middle of the group bars
        plt.xlabel('Material Region', fontproperties=font)
        plt.xticks([r + barWidth/2 for r in range(len(simulated))], labels)

        # Formatting and labels
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
        plt.tick_params(axis='both', which='minor', labelsize=10)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.ylabel('Volume Fraction, $v_a$ (-)', fontproperties=font)
        plt.title('Simulated/Voxelated Material Volume Fraction', fontproperties=font)
        plt.ylim(0,0.1+max(np.concatenate((simulated,voxelated))))

        # Create legend & save graphic
        plt.legend()
        plt.savefig(Path('meshes/' + geoName + '/' + geoName \
            + '-materialVolume.png'), dpi=300)

    def cementMaterialVolPlot(self,geoName,PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,\
        CSH_LDVolFracSim,CSH_HDVolFracSim,PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,\
        CSH_LDVolFracAct,CSH_HDVolFracAct):

        # Plot volume fraction comparison
        plt.close('all')
        fig, ax = plt.subplots()

        labels = ['Pores', 'Clinker', 'CH', 'CSH_LD', 'CSH_HD']

        simulated = [PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim]

        voxelated = [PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct]

        # Set up font
        font = FontProperties()
        font.set_family('serif')
        font.set_name('Times New Roman')
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "dejavuserif"

        # Set width of bar
        barWidth = 0.35
          
        # Set position of bar on X axis
        r1 = np.arange(len(simulated))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
         
        # Make the plot
        plt.bar(r1, simulated, color='0.6', width=barWidth, edgecolor='white', label='Simulated')
        plt.bar(r2, voxelated, color='0.8', width=barWidth, edgecolor='white', label='Voxelated')
        for i, v in enumerate(simulated):
            plt.text(i-0.07, v+0.01, "{:.2f}".format(v), color='k', fontproperties=font)
        for i, v in enumerate(voxelated):
            plt.text(i+0.28, v+0.01, "{:.2f}".format(v), color='k', fontproperties=font)
         
        # Add xticks on the middle of the group bars
        plt.xlabel('Material Region', fontproperties=font)
        plt.xticks([r + barWidth/2 for r in range(len(simulated))], labels)

        # Formatting and labels
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
        plt.tick_params(axis='both', which='minor', labelsize=10)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.ylabel('Volume Fraction, $v_a$ (-)', fontproperties=font)
        plt.title('Simulated/Voxelated Material Volume Fraction', fontproperties=font)
        plt.ylim(0,0.1+max(np.concatenate((simulated,voxelated))))

        # Create legend & save graphic
        plt.legend()
        plt.savefig(Path('meshes/' + geoName + '/' + geoName \
            + '-materialVolume.png'), dpi=300)
        
    def logFile(self,gmshTime,nParticles,placementTime,maxAggD,\
            minAggD,fullerCoef,wcRatio,cementC,volFracAir,q,maxIter,\
            geoName,aggOffset,densityWater,densityCement,allTets,dataType,\
            tetTessTime,writeTime,geoFile,dFiber,lFiber,vFiber,fiberFile,\
            multiMaterial,materialFile,maxGrainD,minGrainD,grainFullerCoef,\
            maxBinderD,minBinderD,binderFullerCoef,maxITZD,minITZD,ITZFullerCoef,output,fibers,\
            itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,binderVolFracAct,\
            aggVolFracAct,sieveCurveDiameter,sieveCurvePassing,matSwitched,materialRule,\
            cementStructure,cementmaterialFile,maxPoresD,minPoresD,PoresFullerCoef,\
            PoresSieveCurveDiameter,PoresSieveCurvePassing,maxClinkerD,minClinkerD,\
            ClinkerFullerCoef,ClinkerSieveCurveDiameter,ClinkerSieveCurvePassing,\
            maxCHD,minCHD,CHFullerCoef,CHSieveCurveDiameter,CHSieveCurvePassing,\
            maxCSH_LDD,minCSH_LDD,CSH_LDFullerCoef,CSH_LDSieveCurveDiameter,CSH_LDSieveCurvePassing,\
            maxCSH_HDD,minCSH_HDD,CSH_HDFullerCoef,CSH_HDSieveCurveDiameter,CSH_HDSieveCurvePassing,\
            PoresVolFracSim,ClinkerVolFracSim,CHVolFracSim,CSH_LDVolFracSim,CSH_HDVolFracSim,\
            PoresVolFracAct,ClinkerVolFracAct,CHVolFracAct,CSH_LDVolFracAct,CSH_HDVolFracAct,outputUnits):

        # Generate log file
        with open(Path('meshes/' + geoName + '/' + geoName + '-report.log'),\
            "w") as f:                                       
            f.write('################################################################################\n')
            f.write('##                Log Created with LDPM Mesh Generation Tool                  ##\n')
            f.write('##             Version 1.6.0 Release - M. Troemner, 13 Mar. 2020              ##\n')
            f.write('################################################################################\n')
            f.write('\n')
            f.write('GEOMETRY:\n') 
            f.write('geoName:                                 ' + geoName + '\n')  
            f.write('\n')
            f.write('CONCRETE PARAMETERS:\n')                   
            f.write('maxAggD:                                 ' + str(maxAggD) + ' ['+ outputUnits + '] \n')  
            f.write('minAggD:                                 ' + str(minAggD) + ' ['+ outputUnits + '] \n')  
            if sieveCurveDiameter == []:
                f.write('fullerCoef:                              ' + str(fullerCoef) + ' [-] \n')
                f.write('q:                                       ' + str(q) + ' [-] \n')  
            else:
                f.write('sieveCurveDiameter:                      ' + str(sieveCurveDiameter) + ' [mm] \n')  
                f.write('sieveCurvePassing:                       ' + str(sieveCurvePassing) + ' [-] \n')  
            f.write('wcRatio:                                 ' + str(wcRatio) + ' [-] \n')  
            f.write('cementC:                                 ' + str(cementC/(1.0E+12)) + ' [tonne/mm3] \n')  
            f.write('volFracAir:                              ' + str(volFracAir) + ' [-] \n')  
            if multiMaterial in ['on','On','Y','y','Yes','yes']: 
                f.write('\nMULTI-MATERIAL PARAMETERS:\n')
                f.write('materialFile:                            ' + str(materialFile) + '\n')
                f.write('materialRule:                            ' + str(materialRule) + '\n')
                if materialRule == 10:
                    f.write('Facets switched:                     ' + str(matSwitched) + '\n')
                f.write('(Aggregate) maxGrainD:                   ' + str(maxGrainD) + ' [mm] \n')
                f.write('(Aggregate) minGrainD:                   ' + str(minGrainD) + ' [mm] \n')
                f.write('(Aggregate) grainFullerCoef:             ' + str(grainFullerCoef) + ' [-] \n')
                f.write('(Aggregate) Simulated Grain Vol Frac:    ' + "{:.5f}".format(aggVolFracSim) + ' [-] \n')
                f.write('(Aggregate) Voxel Grain Vol Frac:        ' + "{:.5f}".format(aggVolFracAct) + ' [-] \n')
                f.write('(Binder)    maxBinderD:                  ' + str(maxBinderD) + ' [mm] \n')
                f.write('(Binder)    minBinderD:                  ' + str(minBinderD) + ' [mm] \n')
                f.write('(Binder)    binderFullerCoef:            ' + str(binderFullerCoef) +' [-] \n')
                f.write('(Binder)    Simulated Binder Vol Frac:   ' + "{:.5f}".format(binderVolFracSim) + ' [-] \n')
                f.write('(Binder)    Voxel Binder Vol Frac:       ' + "{:.5f}".format(binderVolFracAct) + ' [-] \n')
                f.write('(ITZ)       maxITZD:                     ' + str(maxITZD) + ' [mm] \n')
                f.write('(ITZ)       minITZD:                     ' + str(minITZD) + ' [mm] \n')
                f.write('(ITZ)       ITZFullerCoef:               ' + str(ITZFullerCoef) + ' [-] \n')
                f.write('(ITZ)       Simulated ITZ Vol Frac:      ' + "{:.5f}".format(itzVolFracSim) + ' [-] \n')
                f.write('(ITZ)       Voxel ITZ Vol Frac:          ' + "{:.5f}".format(itzVolFracAct) + ' [-] \n')

            if cementStructure in ['on','On','Y','y','Yes','yes']: 
                f.write('\nCEMENTSTRUCTURE PARAMETERS:\n')
                f.write('cementmaterialFile:                      ' + str(cementmaterialFile) + '\n')
                f.write('materialRule:                            ' + str(materialRule) + '\n')
                if materialRule == 10:
                    f.write('Facets switched:                     ' + str(matSwitched) + '\n')
                f.write('(Pores) maxPoresD:                       ' + str(maxPoresD) + ' [mm] \n')
                f.write('(Pores) minPoresD:                       ' + str(minPoresD) + ' [mm] \n')
                f.write('(Pores) PoresFullerCoef:                 ' + str(PoresFullerCoef) + ' [-] \n')
                f.write('(Pores) Simulated Pores Vol Frac:        ' + "{:.5f}".format(PoresVolFracSim) + ' [-] \n')
                f.write('(Pores) Voxel Pores Vol Frac:            ' + "{:.5f}".format(PoresVolFracAct) + ' [-] \n')
                f.write('(Clinker)   maxClinkerD:                 ' + str(maxClinkerD) + ' [mm] \n')
                f.write('(Clinker)   minClinkerD:                 ' + str(minClinkerD) + ' [mm] \n')
                f.write('(Clinker)   ClinkerFullerCoef:           ' + str(ClinkerFullerCoef) +' [-] \n')
                f.write('(Clinker)   Simulated Clinker Vol Frac   ' + "{:.5f}".format(ClinkerVolFracSim) + ' [-] \n')
                f.write('(Clinker)   Voxel Clinker Vol Frac:      ' + "{:.5f}".format(ClinkerVolFracAct) + ' [-] \n')
                f.write('(CH)       maxCHD:                       ' + str(maxCHD) + ' [mm] \n')
                f.write('(CH)       minCHD:                       ' + str(minCHD) + ' [mm] \n')
                f.write('(CH)       CHFullerCoef:                 ' + str(CHFullerCoef) + ' [-] \n')
                f.write('(CH)       Simulated CH Vol Frac:        ' + "{:.5f}".format(CHVolFracSim) + ' [-] \n')
                f.write('(CH)       Voxel CH Vol Frac:            ' + "{:.5f}".format(CHVolFracAct) + ' [-] \n')
                f.write('(CSH_LD)  maxCSH_LDD:                    ' + str(maxCSH_LDD) + ' [mm] \n')
                f.write('(CSH_LD)  minCSH_LDD:                    ' + str(minCSH_LDD) + ' [mm] \n')
                f.write('(CSH_LD)  CSH_LDFullerCoef:              ' + str(CSH_LDFullerCoef) + ' [-] \n')
                f.write('(CSH_LD)  Simulated CSH_LD Vol Frac:     ' + "{:.5f}".format(CSH_LDVolFracSim) + ' [-] \n')
                f.write('(CSH_LD)  Voxel CSH_LD Vol Frac:         ' + "{:.5f}".format(CSH_LDVolFracAct) + ' [-] \n')
                f.write('(CSH_HD)  maxCSH_HDD:                    ' + str(maxCSH_HDD) + ' [mm] \n')
                f.write('(CSH_HD)  minCSH_HDD:                    ' + str(minCSH_HDD) + ' [mm] \n')
                f.write('(CSH_HD)  CSH_HDFullerCoef:              ' + str(CSH_HDFullerCoef) + ' [-] \n')
                f.write('(CSH_HD)  Simulated CSH_HD Vol Frac:     ' + "{:.5f}".format(CSH_HDVolFracSim) + ' [-] \n')
                f.write('(CSH_HD)  Voxel CSH_HD Vol Frac:         ' + "{:.5f}".format(CSH_HDVolFracAct) + ' [-] \n')

            f.write('\nGENERATION PARAMETERS:\n')
            f.write('maxIter:                                 ' + str(maxIter) + ' [-] \n')  
            f.write('dataType:                                ' + dataType + '\n')
            f.write('\n')
            f.write('GENERATION:\n')  
            f.write('Number of Particles:                     ' + str(nParticles) + '\n') 
            f.write('Number of Tets:                          ' + str(len(allTets)) + '\n') 
            f.write('Number of Facets:                        ' + str(len(allTets)*12) + '\n') 
            f.write('\n')
            f.write('PERFORMANCE:\n')  
            f.write('Gmsh time:                               ' + str(gmshTime) + ' [seconds] \n')  
            f.write('Placement Time:                          ' + str(placementTime) + ' [seconds] \n')     
            f.write('Tetrahedralization/Tesselation time:     ' + str(tetTessTime) \
                + ' [seconds] \n')  
            f.write('Writing Time:                            ' + str(writeTime) + ' [seconds] \n')     
            f.write('Total Time:                              ' + str(round(time.time()-start_time,2)) \
                + ' [seconds] \n')          
            f.write('\n')
            f.write('FILE CREATION:\n')  
            f.write('Files Requested:                         ' + str(output) + '\n')
            f.write('Geometry File:                           ./meshes/' + geoName + '/' + geoFile + '\n')
            f.write('Cell Data File:                          ./meshes/' + geoName + '/' + geoName \
                + '-data-cell.dat\n')
            f.write('Facet Data File:                         ./meshes/' + geoName + '/' + geoName \
                + '-data-facet.dat\n')
            f.write('Mesh Data File:                          ./meshes/' + geoName + '/' + geoName \
                + '-data-mesh.dat\n')
            f.write('Paraview Mesh File:                      ./meshes/' + geoName + '/' + geoName
                + '-para-mesh.000.vtk\n')
            f.write('Paraview Particles File:                 ./meshes/' + geoName + '/' \
                + geoName + '-para-particles.000.vtk\n')
            f.write('Paraview Facet File:                     ./meshes/' + geoName + '/' + geoName \
                + '-para-facet.000.vtk\n')
            f.write('Log File:                                ./meshes/' + geoName + '/' + geoName \
                + '-report.log\n')
            if fibers in ['on','On','Y','y','Yes','yes']:
                f.write('Fiber Facet Data File:                   ./meshes/' + geoName + '/' + geoName \
                    + '-data-fiberfacet.dat\n')

    # Section only used for study and debugging; leave as-is commented out.

    #def materialLogFile(self,materialFile,itzVolFracSim,binderVolFracSim,aggVolFracSim,itzVolFracAct,\
    #   binderVolFracAct,aggVolFracAct,materialRule,minAggD):

        # Generate material fraction log file
    #   with open(Path('meshes/' + materialFile + '-MaterialReport.log'),\
    #       "a") as f:                                       
    #       f.write('Rule: ' + str(materialRule) + ' minAggD: ' + str(minAggD) + ' aggVolFracSim: ' + "{:.8f}".format(aggVolFracSim) + ' aggVolFracAct: ' + "{:.8f}".format(aggVolFracAct) + ' binderVolFracSim: ' + "{:.8f}".format(binderVolFracSim) + ' binderVolFracAct: ' + "{:.8f}".format(binderVolFracAct) + ' itzVolFracSim: ' + "{:.8f}".format(itzVolFracSim) + ' itzVolFracAct: ' + "{:.8f}".format(itzVolFracAct) + '\n')
