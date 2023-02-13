from pathlib import Path


def mkDataNodes(geoName,tempPath,allNodes):
    
    with open(Path(tempPath + geoName + \
        '-data-nodes.dat'),"w") as f:                                            
        f.write('// Mesh Data Generated with LDPM Mesh Generation Tool\n')
        f.write('\n') 
        f.write('// Number of Vertices: ' + str(len(allNodes)) + '\n')  
        f.write("\n".join(" ".join(map(str, x)) for x in allNodes))     
        f.write('\n')