## ===========================================================================
## CHRONO WORKBENCH:github.com/Concrete-Chrono-Development/chrono-preprocessor
##
## Copyright (c) 2023 
## All rights reserved. 
##
## Use of this source code is governed by a BSD-style license that can be
## found in the LICENSE file at the top level of the distribution and at
## github.com/Concrete-Chrono-Development/chrono-preprocessor/blob/main/LICENSE
##
## ===========================================================================
## Developed by Northwestern University
## For U.S. Army ERDC Contract No. W9132T22C0015
## Primary Authors: Matthew Troemner
## ===========================================================================
##
## This function creates a 3D geometric shape using dimensions and parameters 
## passed to it. The function supports several different shapes including 
## boxes, cylinders, cones, spheres, ellipsoids, arbitrary prisms, and notched 
## prisms of square, semi-circle, and semi-ellipse shapes. The 'Dogbone' 
## option creates a special 3D dogbone shape using specific dimensions passed 
## to the function.
##
## ===========================================================================

# pyright: reportMissingImports=false
import os
import re
import math

import FreeCAD as App
import ImportGui
import JoinFeatures
import BOPTools.JoinFeatures
import Part
from FreeCAD import Base


def gen_LDPMCSL_geometry(dimensions,geoType,geoName,cadFile):

    """
    Variables:
    --------------------------------------------------------------------------
    ### Inputs ###
    - dimensions: List of dimensions for the geometry
    - geoType: Type of geometry to be created
    - geoName: Name of the geometry
    - cadFile: Path to the CAD file to be imported
    --------------------------------------------------------------------------
    ### Outputs ###
    - geo: Geometry object
    --------------------------------------------------------------------------
    """  


    # Check if dimensions are positive (ignore geometries we cannot check)
    if geoType not in ['Ellipsoid', 'Dogbone', 'Custom', 'Import CAD']:
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






    if geoType == "Dogbone":

        # Define variables
        length = float(dimensions[0].strip(" mm"))  # Length of the part
        width = float(dimensions[1].strip(" mm"))  # Width of the part
        thickness = float(dimensions[2].strip(" mm"))  # Thickness of the part
        gauge_length = float(dimensions[3].strip(" mm"))  # Gauge length of the dogbone
        gauge_width = float(dimensions[4].strip(" mm"))  # Gauge width of the dogbone
        dogbone_type = dimensions[5]  # Type of the dogbone shape

        # Create a rectangular shape
        rectangle = App.ActiveDocument.addObject("Part::Box", geoName+"Rectangle")
        rectangle.Length = width
        rectangle.Width = thickness
        rectangle.Height = length

        # Create a dogbone shape
        if dogbone_type == 'Rounded':
            radius = (width-gauge_width)/2
            dogbone1 = App.ActiveDocument.addObject("Part::Cylinder", geoName+"Dogbone1")
            dogbone1.Radius = radius
            dogbone1.Height = thickness
            dogbone1.Placement = App.Placement(App.Vector(0, 0, length/2-gauge_length/2),App.Rotation(App.Vector(1,0,0),-90.00))

            dogbone2 = App.ActiveDocument.addObject("Part::Cylinder", geoName+"Dogbone2")
            dogbone2.Radius = radius
            dogbone2.Height = thickness
            dogbone2.Placement = App.Placement(App.Vector(width, 0, length/2-gauge_length/2),App.Rotation(App.Vector(1,0,0),-90.00))

            dogbone3 = App.ActiveDocument.addObject("Part::Cylinder", geoName+"Dogbone3")
            dogbone3.Radius = radius
            dogbone3.Height = thickness
            dogbone3.Placement = App.Placement(App.Vector(0, 0, length/2+gauge_length/2),App.Rotation(App.Vector(1,0,0),-90.00))

            dogbone4 = App.ActiveDocument.addObject("Part::Cylinder", geoName+"Dogbone4")
            dogbone4.Radius = radius
            dogbone4.Height = thickness
            dogbone4.Placement = App.Placement(App.Vector(width, 0, length/2+gauge_length/2),App.Rotation(App.Vector(1,0,0),-90.00))

            side1 = App.ActiveDocument.addObject("Part::Box", geoName+"side1")
            side1.Length = (width-gauge_width)/2
            side1.Width = thickness
            side1.Height = gauge_length
            side1.Placement = App.Placement(App.Vector(0, 0, length/2-gauge_length/2),App.Rotation(App.Vector(0,0,0),0.00))

            side2 = App.ActiveDocument.addObject("Part::Box", geoName+"side2")
            side2.Length = (width-gauge_width)/2
            side2.Width = thickness
            side2.Height = gauge_length
            side2.Placement = App.Placement(App.Vector(width-(width-gauge_width)/2, 0, length/2-gauge_length/2),App.Rotation(App.Vector(0,0,0),0.00))

            DogboneInsert = App.ActiveDocument.addObject("Part::MultiFuse", geoName+"DogboneInsert")
            DogboneInsert.Shapes = [dogbone1, dogbone2, dogbone3, dogbone4, side1, side2]

        else:
            side1 = App.ActiveDocument.addObject("Part::Box", geoName+"side1")
            side1.Length = (width-gauge_width)/2
            side1.Width = thickness
            side1.Height = gauge_length
            side1.Placement = App.Placement(App.Vector(0, 0, length/2-gauge_length/2),App.Rotation(App.Vector(0,0,0),0.00))

            side2 = App.ActiveDocument.addObject("Part::Box", geoName+"side2")
            side2.Length = (width-gauge_width)/2
            side2.Width = thickness
            side2.Height = gauge_length
            side2.Placement = App.Placement(App.Vector(width-(width-gauge_width)/2, 0, length/2-gauge_length/2),App.Rotation(App.Vector(0,0,0),0.00))

            DogboneInsert = App.ActiveDocument.addObject("Part::MultiFuse", geoName+"DogboneInsert")
            DogboneInsert.Shapes = [side1, side2]

        # Subtract the dogbone from the rectangle
        part = App.ActiveDocument.addObject("Part::Cut",geoName)
        part.Base = rectangle
        part.Tool = DogboneInsert


    if geoType == "Custom":
        # Get part name from self.form[1].selectedObject
        #partName = form[1].selectedObject.toPlainText()

        # Find the part in the document and change its name to geoName
        #part = App.getDocument(App.ActiveDocument.Name).getObjectsByLabel(partName)[0]
        #part.Label = geoName
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
