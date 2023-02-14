import os
import re

import FreeCAD as App
import ImportGui
import JoinFeatures
import BOPTools.JoinFeatures
from FreeCAD import Base


def genGeometry(dimensions,geoType,geoName,cadFile):

    # Check if dimensions are positive (ignore geometries we cannot check)
    if geoType not in ['Ellipsoid', 'Custom', 'Import CAD']:
        if all(float(i.strip(" mm")) > 0 for i in dimensions):
            pass
        else:
            raise Exception("One or more geometry dimensions are less than or equal to zero. Please revise.")




    if geoType == "Box":

        # Create a box and name it
        geo           = App.ActiveDocument.addObject("Part::Box",geoName)
        geo.Label     = geoName
        geo.Height    = dimensions[0]
        geo.Width     = dimensions[1]
        geo.Length    = dimensions[2]


    if geoType == "Cylinder":

        # Create a box and name it
        geo           = App.ActiveDocument.addObject("Part::Cylinder",geoName)
        geo.Label     = geoName
        geo.Height    = dimensions[0]
        geo.Radius    = dimensions[1]

    if geoType == "Cone":

        # Create a cone and name it
        geo           = App.ActiveDocument.addObject("Part::Cone",geoName)
        geo.Label     = geoName
        geo.Height    = dimensions[0]
        geo.Radius1   = dimensions[1]
        geo.Radius2   = dimensions[2]

    if geoType == "Sphere":

        # Create a sphere and name it
        geo           = App.ActiveDocument.addObject("Part::Sphere",geoName)
        geo.Label     = geoName
        geo.Radius    = dimensions[0]

    if geoType == "Ellipsoid":

        # Create an ellipsoid and name it
        geo           = App.ActiveDocument.addObject("Part::Ellipsoid",geoName)
        geo.Label     = geoName
        geo.Radius1   = dimensions[0]
        geo.Radius2   = dimensions[1]
        geo.Radius3   = dimensions[2]
        geo.Angle1    = dimensions[3]
        geo.Angle2    = dimensions[4]
        geo.Angle3    = dimensions[5]

    if geoType == "Arbitrary Prism":

        # Create a prism and name it
        geo           = App.ActiveDocument.addObject("Part::Prism",geoName)
        geo.Label     = geoName
        geo.Circumradius = dimensions[0]
        geo.Height    = dimensions[1]
        geo.Polygon   = int(dimensions[2])

    if geoType == "Notched Prism - Square":

        # Create a notched prism and name it
        geo           = App.ActiveDocument.addObject("Part::Box",geoName+"Box")
        geo.Length    = dimensions[0]
        geo.Width     = dimensions[1]
        geo.Height    = dimensions[2]
        geo.Placement = App.Placement(App.Vector(0.00,0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))
        geo.Label     = geoName+'Box'

        # Create the notch
        geo           = App.ActiveDocument.addObject("Part::Box",geoName+"Notch")
        geo.Length    = dimensions[3]    # Notch Width
        geo.Width     = dimensions[1]    # Box width
        geo.Height    = dimensions[4]    # Notch Depth
        geo.Placement = App.Placement(App.Vector(0.00, 0.00, 0.00), App.Rotation(App.Vector(0.00, 0.00, 1.00), 0.00))
        geo.Label     = geoName + 'Notch'
        App.getDocument(App.ActiveDocument.Name).getObject(geoName+'Notch').Placement = App.Placement(App.Vector((float(dimensions[0].strip(" mm"))/2-float(dimensions[3].strip(" mm"))/2),0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))

        # Cut out the notch
        j = BOPTools.JoinFeatures.makeCutout(name=geoName)
        j.Base = App.ActiveDocument.getObject(geoName+'Box')
        j.Tool = App.ActiveDocument.getObject(geoName+'Notch')
        j.Proxy.execute(j)
        j.purgeTouched()
        for obj in j.ViewObject.Proxy.claimChildren():
            obj.ViewObject.hide()

    if geoType == "Notched Prism - Semi Circle":

        # Create a notched prism and name it
        geo            = App.ActiveDocument.addObject("Part::Box",geoName+"Box")
        geo.Length     = dimensions[0]
        geo.Width      = dimensions[1]
        geo.Height     = dimensions[2]
        geo.Placement  = App.Placement(App.Vector(0.00,0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))
        geo.Label      = geoName+'Box'

        geo            = App.ActiveDocument.addObject("Part::Box",geoName+"BoxNotch")
        geo.Length     = dimensions[3] # BoxNotch Width
        geo.Width      = dimensions[1]
        geo.Height     = dimensions[4] # BoxNotch Depth
        geo.Placement  = App.Placement(App.Vector(float(dimensions[0].strip(" mm"))/2-float(dimensions[3].strip(" mm"))/2,0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))
        geo.Label      = geoName+'BoxNotch'

        geo            = App.ActiveDocument.addObject("Part::Cylinder",geoName+"CylinderNotch")
        geo.Radius     = float(dimensions[3].strip(" mm"))/2
        geo.Height     = float(dimensions[1].strip(" mm"))
        geo.Angle      = 360.00
        geo.FirstAngle = 0.00
        geo.SecondAngle= 0.00
        geo.Placement  = App.Placement(App.Vector(float(dimensions[0].strip(" mm"))/2,0.00,float(dimensions[4].strip(" mm"))),App.Rotation(App.Vector(1.00,0.00,0.00),-90.00))
        geo.Label      = geoName+'CylinderNotch'

        j = BOPTools.JoinFeatures.makeConnect(name=geoName+'Connect')
        j.Objects = [App.ActiveDocument.getObject(geoName+'BoxNotch'), App.ActiveDocument.getObject(geoName+'CylinderNotch')]
        j.Proxy.execute(j)
        j.purgeTouched()
        for obj in j.ViewObject.Proxy.claimChildren():
            obj.ViewObject.hide()

        j = BOPTools.JoinFeatures.makeCutout(name=geoName)
        j.Base = App.ActiveDocument.getObject(geoName+'Box')
        j.Tool = App.ActiveDocument.getObject(geoName+'Connect')
        j.Proxy.execute(j)
        j.purgeTouched()
        for obj in j.ViewObject.Proxy.claimChildren():
            obj.ViewObject.hide()



    if geoType == "Notched Prism - Semi Ellipse":

        # Create a notched prism and name it
        geo            = App.ActiveDocument.addObject("Part::Box",geoName+"Box")
        geo.Length     = dimensions[0]
        geo.Width      = dimensions[1]
        geo.Height     = dimensions[2]
        geo.Placement  = App.Placement(App.Vector(0.00,0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))
        geo.Label      = geoName+'Box'

        geo            = App.ActiveDocument.addObject("Part::Box",geoName+"BoxNotch")
        geo.Length     = dimensions[3] # BoxNotch Width
        geo.Width      = dimensions[1]
        geo.Height     = dimensions[4] # BoxNotch Depth
        geo.Placement  = App.Placement(App.Vector(float(dimensions[0].strip(" mm"))/2-float(dimensions[3].strip(" mm"))/2,0.00,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))
        geo.Label      = geoName+'BoxNotch'

        if float(dimensions[5].strip(" mm")) > float(dimensions[3].strip(" mm"))/2:
            geo        = App.ActiveDocument.addObject("Part::Ellipse",geoName+"Ellipse")
            geo.MajorRadius = float(dimensions[5].strip(" mm"))
            geo.MinorRadius = float(dimensions[3].strip(" mm"))/2
            geo.Angle1 = 0.00
            geo.Angle2 = 360.00
            geo.Placement = App.Placement(App.Vector(0.00,0.00,0.00),App.Rotation(App.Vector(1.00,0.00,0.00),-90.00))
            geo.Label  = geoName+"Ellipse"
        else:
            geo        = App.ActiveDocument.addObject("Part::Ellipse",geoName+"Ellipse")
            geo.MajorRadius = float(dimensions[3].strip(" mm"))/2
            geo.MinorRadius = float(dimensions[5].strip(" mm"))
            geo.Angle1 = 0.00
            geo.Angle2 = 360.00
            geo.Placement = App.Placement(App.Vector(0.00,0.00,0.00),App.Rotation(App.Vector(1.00,0.00,0.00),-90.00))
            geo.Label  = geoName+"Ellipse"

        f              = App.ActiveDocument.addObject('Part::Extrusion',geoName+'EllipseNotch')
        f.Base         = App.ActiveDocument.getObject(geoName+'Ellipse')
        f.DirMode      = "Normal"
        f.DirLink      = None
        f.LengthFwd    = float(dimensions[1].strip(" mm"))
        f.LengthRev    = 0.000000000000000
        f.Solid        = True
        f.Reversed     = False
        f.Symmetric    = False
        f.TaperAngle   = 0.000000000000000
        f.TaperAngleRev= 0.000000000000000

        if float(dimensions[5].strip(" mm")) > float(dimensions[3].strip(" mm"))/2:
            App.ActiveDocument.getObject(geoName+'EllipseNotch').Placement = App.Placement(App.Vector(float(dimensions[0].strip(" mm"))/2,0.00,float(dimensions[4].strip(" mm"))),App.Rotation(App.Vector(0.00,1.00,0.00),90.00))
        else:
            App.ActiveDocument.getObject(geoName+'EllipseNotch').Placement = App.Placement(App.Vector(float(dimensions[0].strip(" mm"))/2,0.00,float(dimensions[4].strip(" mm"))),App.Rotation(App.Vector(0.00,1.00,0.00),0.00))

        App.ActiveDocument.recompute()

        j = BOPTools.JoinFeatures.makeConnect(name=geoName+'Connect')
        j.Objects = [App.ActiveDocument.getObject(geoName+'BoxNotch'), App.ActiveDocument.getObject(geoName+'EllipseNotch')]
        j.Proxy.execute(j)
        j.purgeTouched()
        for obj in j.ViewObject.Proxy.claimChildren():
            obj.ViewObject.hide()

        j = BOPTools.JoinFeatures.makeCutout(name=geoName)
        j.Base = App.ActiveDocument.getObject(geoName+'Box')
        j.Tool = App.ActiveDocument.getObject(geoName+'Connect')
        j.Proxy.execute(j)
        j.purgeTouched()
        for obj in j.ViewObject.Proxy.claimChildren():
            obj.ViewObject.hide()

        App.ActiveDocument.getObject(geoName+'Ellipse').Visibility = False
        App.ActiveDocument.getObject(geoName+'EllipseNotch').Visibility = False
        App.ActiveDocument.getObject(geoName+'BoxNotch').Visibility = False



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