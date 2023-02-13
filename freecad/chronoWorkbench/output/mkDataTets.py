from pathlib import Path


def mkDataTets(geoName,tempPath,allTets):
    
    with open(Path(tempPath + geoName + \
        '-data-tets.dat'),"w") as f:                                            
        f.write('// Mesh Data Generated with LDPM Mesh Generation Tool\n')
        f.write('\n') 
        f.write('// Number of Tets: ' + str(len(allTets)) + '\n')  
        f.write("\n".join(" ".join(map(str, x)) for x in \
            allTets.astype(int)))   
        f.write('\n')