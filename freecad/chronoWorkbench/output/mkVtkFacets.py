import numpy as np
from pathlib import Path





def mkVtkFacets(geoName,tetFacets,dataList,facetMaterial,\
    multiMaterial,cementStructure,tempPath,facetPointData,facetCellData):

    #facetsPoints = tetFacets.reshape(-1,3)
    #cells = (np.arange(0,round(len(facetsPoints))).\
    #    reshape(-1,3)).astype(int)
    #cell_types = np.array([5,]*round(len(facetsPoints)/3))
    cell_types = np.array([5,]*(len(facetCellData)))
    while facetPointData.shape[0] % 3 != 0:
        facetPointData = np.concatenate((facetPointData, np.zeros((1, facetPointData.shape[1]))), axis=0)

    facetPointData = np.around(facetPointData.reshape(-1,9), decimals=6) # reshape and condense to save memory/space

    with open(Path(tempPath + geoName + \
        '-para-facet.000.vtk'),"w") as f:                                                                          
        f.write('# vtk DataFile Version 2.0\n')
        f.write('Facet Visual File\n') 
        f.write('ASCII\n')    
        f.write('DATASET POLYDATA\n')        


        #f.write('POINTS ' + str(len(facetsPoints)) + ' double \n')  
        #f.write("\n".join(" ".join(map(str, x)) for x in facetsPoints))
        f.write('POINTS ' + str(len(facetPointData)*3) + ' float \n') 
        f.write("\n".join(" ".join(map(str, x)) for x in facetPointData))
        f.write('\n\n')  
        #f.write('CELLS ' + str(round(len(facetsPoints)/3)) + ' ' \
        #    + str(round(len(facetsPoints)/3*4)) +'\n3 ')
        #f.write("\n3 ".join(" ".join(map(str, x)) for x in cells))
        f.write('POLYGONS ' + str(len(facetCellData)) + ' ' \
            + str(round(len(facetCellData)*4)) +'\n3 ')
        f.write("\n3 ".join(" ".join(map(str, x)) for x in facetCellData))

        f.write('\n\n')  
        #f.write('CELL_TYPES ' + str(round(len(facetCellData))) +'\n')
        #for x in cell_types:
        #    f.write("%s\n" % x)
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