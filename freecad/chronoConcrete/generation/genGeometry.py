
import FreeCAD as App


def genGeometry(dimensions,geoType,geoName):


    if all(float(i.strip(" mm")) > 0 for i in dimensions):
        pass
    else:
        raise Exception("One or more geometry dimensions are less than or equal to zero. Please revise.")




    if geoType == "Box":

        # Create a box and name it
        geo = App.ActiveDocument.addObject("Part::Box",geoName)
        geo.Label = geoName
        geo.Height = dimensions[0]
        geo.Width = dimensions[1]
        geo.Length = dimensions[2]


    App.ActiveDocument.recompute()
