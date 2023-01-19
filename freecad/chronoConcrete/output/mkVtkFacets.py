import numpy as np
from pathlib import Path





def mkVtkFacets(geoName,tetFacets,dataList,facetMaterial,\
    multiMaterial,cementStructure,tempPath):

    facetsPoints = tetFacets.reshape(-1,3)
    cells = (np.arange(0,round(len(facetsPoints))).\
        reshape(-1,3)).astype(int)
    cell_types = np.array([5,]*round(len(facetsPoints)/3))

    with open(Path(tempPath + geoName + \
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