
import vtk

def extract_edges_from_vtk(input_file, output_file):
    # Read the input VTK file
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(input_file)
    reader.Update()
    
    # Extract the edges from the input data
    edges = vtk.vtkExtractEdges()
    edges.SetInputData(reader.GetOutput())
    edges.Update()
    
    # Write the output VTK file
    writer = vtk.vtkDataSetWriter()
    writer.SetFileVersion(4)
    writer.SetFileName(output_file)
    writer.SetFileTypeToASCII()
    writer.SetInputData(edges.GetOutput())
    writer.Write()
