import FreeCAD as App



# The `genGeometry` function creates a geometric object of specified type and dimensions and adds it to the current FreeCAD document.

def genGeometry(dimensions, geoType, geoName):
    # Check if the given geometry type is valid
    if geoType not in ["Box", "Cylinder", "Cone", "Sphere", "Ellipsoid", "Prism"]:
        raise Exception("Invalid geometry type. Please choose from 'Box', 'Cylinder', 'Cone', 'Sphere', 'Ellipsoid', 'Prism'.")

    # Check if all dimensions are positive
    if any(float(i.strip(" mm")) <= 0 for i in dimensions):
        raise Exception("One or more geometry dimensions are less than or equal to zero. Please revise.")

    # Create the specified geometry object
    if geoType == "Box":
        geo = App.ActiveDocument.addObject("Part::Box", geoName) # Add a box to the document
        geo.Label = geoName # Set the label of the box
        geo.Height, geo.Width, geo.Length = dimensions # Set the dimensions of the box
    elif geoType == "Cylinder":
        geo = App.ActiveDocument.addObject("Part::Cylinder", geoName) # Add a cylinder to the document
        geo.Label = geoName # Set the label of the cylinder
        geo.Height, geo.Radius = dimensions # Set the dimensions of the cylinder
    elif geoType == "Cone":
        geo = App.ActiveDocument.addObject("Part::Cone", geoName) # Add a cone to the document
        geo.Label = geoName # Set the label of the cone
        geo.Height, geo.Radius1, geo.Radius2 = dimensions # Set the dimensions of the cone
    elif geoType == "Sphere":
        geo = App.ActiveDocument.addObject("Part::Sphere", geoName) # Add a sphere to the document
        geo.Label = geoName # Set the label of the sphere
        geo.Radius = dimensions[0] # Set the radius of the sphere
    elif geoType == "Ellipsoid":
        geo = App.ActiveDocument.addObject("Part::Ellipsoid", geoName) # Add an ellipsoid to the document
        geo.Label = geoName # Set the label of the ellipsoid
        geo.Radius1, geo.Radius2, geo.Radius3, geo.Angle1, geo.Angle2, geo.Angle3 = dimensions # Set the dimensions of the ellipsoid
    elif geoType == "Prism":
        geo = App.ActiveDocument.addObject("Part::Prism", geoName) # Add a prism to the document
        geo.Label = geoName # Set the label of the prism
        geo.Circumradius, geo.Height, geo.Polygon = dimensions # Set the dimensions of the prism
        geo.Polygon = int(geo.Polygon) # Convert the polygon sides to an integer

    # Recompute the document to reflect the changes
    App.ActiveDocument.recompute()
