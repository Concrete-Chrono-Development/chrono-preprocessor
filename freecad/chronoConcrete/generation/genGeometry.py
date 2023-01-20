
import FreeCAD as App


def genGeometry(dimensions,geoType,geoName):


    if geoType in ['Box', 'Cylinder', 'Cone', 'Sphere']:
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


    if geoType == "Cylinder":

        # Create a box and name it
        geo = App.ActiveDocument.addObject("Part::Cylinder",geoName)
        geo.Label = geoName
        geo.Height = dimensions[0]
        geo.Radius = dimensions[1]

    if geoType == "Cone":

        # Create a cone and name it
        geo = App.ActiveDocument.addObject("Part::Cone",geoName)
        geo.Label = geoName
        geo.Height = dimensions[0]
        geo.Radius1 = dimensions[1]
        geo.Radius2 = dimensions[2]

    if geoType == "Sphere":

        # Create a sphere and name it
        geo = App.ActiveDocument.addObject("Part::Sphere",geoName)
        geo.Label = geoName
        geo.Radius = dimensions[0]

    if geoType == "Ellipsoid":

        # Create an ellipsoid and name it
        geo = App.ActiveDocument.addObject("Part::Ellipsoid",geoName)
        geo.Label = geoName
        geo.Radius1 = dimensions[0]
        geo.Radius2 = dimensions[1]
        geo.Radius3 = dimensions[2]
        geo.Angle1 = dimensions[3]
        geo.Angle2 = dimensions[4]
        geo.Angle3 = dimensions[5]

    if geoType == "Prism":

        # Create a prism and name it
        geo = App.ActiveDocument.addObject("Part::Prism",geoName)
        geo.Label = geoName
        geo.Circumradius = dimensions[0]
        geo.Height = dimensions[1]
        geo.Polygon = int(dimensions[2])




    App.ActiveDocument.recompute()
