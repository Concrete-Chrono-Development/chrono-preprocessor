import os
import re

import FreeCAD as App
import ImportGui
import JoinFeatures
import BOPTools.JoinFeatures


def genGeometry(dimensions,geoType,geoName,cadFile):


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

    if geoType == "Notched Prism":

        # Create a notched prism and name it
        App.ActiveDocument.addObject("Part::Box","Box")
        App.ActiveDocument.Box.Length=dimensions[0]
        App.ActiveDocument.Box.Width=dimensions[1]
        App.ActiveDocument.Box.Height=dimensions[2]
        App.ActiveDocument.Box.Placement=App.Placement(App.Vector(0.00,0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))
        App.ActiveDocument.Box.Label='Box'

        App.ActiveDocument.addObject("Part::Box","Notch")
        App.ActiveDocument.Notch.Length=dimensions[3] # Notch Width
        App.ActiveDocument.Notch.Width=dimensions[1]
        App.ActiveDocument.Notch.Height=dimensions[4] # Notch Depth
        App.ActiveDocument.Notch.Placement=App.Placement(App.Vector(0.00,0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))
        App.ActiveDocument.Notch.Label='Notch'

        App.getDocument(App.ActiveDocument.Name).getObject('Notch').Placement = App.Placement(App.Vector((float(dimensions[0].strip(" mm"))/2-float(dimensions[3].strip(" mm"))/2),0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))

        j = BOPTools.JoinFeatures.makeCutout(name=geoName)
        j.Base = App.ActiveDocument.Box
        j.Tool = App.ActiveDocument.Notch
        j.Proxy.execute(j)
        j.purgeTouched()
        for obj in j.ViewObject.Proxy.claimChildren():
            obj.ViewObject.hide()


    if geoType == "Custom":
        pass


    if geoType == "Import CAD":
        ImportGui.insert(cadFile,App.ActiveDocument.Name)
        filename = os.path.basename(cadFile)
        filename, file_extension = os.path.splitext(filename)
        filename = re.sub("\.", "_", filename)
        filename = re.sub("/.", "_", filename)
        geo = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(filename)[0]
        geo.Label = geoName

    App.ActiveDocument.recompute()
