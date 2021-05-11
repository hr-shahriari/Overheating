
"""Use this plugin to take any list of closed brep masses to calculate the risk of overheating. By Hamidreza Shahriari
    Inputs:
        North: A vector to represent the North direction. If not provided the default Y axis will be selected as North
        Latitude: Latitude of the locatation of the building
        masses: A list of closed Breps to perform the overheating study on the building masses
        g_value: g_value of the windows
        bld_type: Building type according to the SAP Table p1. The values: {"single story dwelling with possibility of cross ventilation" = 0 , "single story dwelling cross ventilation not possible" = 1, "Dwelling of two or more storeys windows open upstairs and downstairs Cross ventilation possible" = 2, "Dwelling of two or more storeys windows open upstairs and downstairs Cross ventilation not possible" =3} 
        solar_access: Percentage of Sky blocked by obstecles. {"Heavy" > 80, "More than avarage" >60 & <80, "average or Unknown" >20 & <60, "very little" <20} 
        U_window: Total U_value of the window system
        U_wall : U_value of the walls.
        WWR: Window to WAll Ratio
        solar_flux: solar flux on the horizontal surface during the summer period
        Roof: Set to True if there is roof
        max_roof_angle: Maximum angle of the pitched roof
        U_roof: U_value of the roof
        overhangs: Overhang depth from the glass (in meter)
        curt_blind: Selcet the type of blinds, curtains or external shutters to calculate the shading factor according to SAP table P3: {"Net curtain (covering whole window)" = 0, " Net curtain (covering half window)" = 1, " Dark-coloured curtain or roller blind (note 1)" = 2, " Light-coloured curtain or roller blind (note 1)" = 3, " Dark-coloured venetian blind (note 2)" = 4," Light-coloured venetian blind (note 2)” = 5, “Dark-coloured external shutter, window closed (notes 1, 3)" = 6, “White external shutter, window closed (notes 1, 3)” = 7, “Dark-coloured external shutter, window fully open (notes 1, 3)” = 8, “White external shutter, window fully open (notes 1, 3)” = 9}
        opn_type: Ventilation Opening type according to the SAP Table p1.he values: {"Trickle vents only" = 0 , "Windows slightly open (50 mm)" = 1, "Windows open half the time" = 2,"Windows fully open" = 3}
        TMP: Thermal mass parameter
        FF: The frame factor for windows and doors (fraction of opening that is glazed). Frame factor according to material: {"wood" = 0.7, "metal" = 0.8, PVC-U = 0.7} 
        ex_t_average: The mean external temperature for the selected summer month and climate region.
        time_fraction: fraction of the daylight hours that shading blinds are in use
        window_height: The height of the window frame
        month_Solar_declination: 
    Outputs:
        results: Likelihood of high internal temperature during hot weather.
        preview: Colored mesh preview
"""
#Copyright (c) 2021, Hamidreza Shahriari 
__author__ = "Hamidreza"



import rhinoscriptsyntax as rs
import Rhino.Geometry as rg
import ghpythonlib.components as ghc
import math
import Grasshopper.Kernel.GH_Convert as ghconvert
from scriptcontext import doc

if not max_roof_angle:
    max_roof_angle = 1

def air_change_rate(op_type,bld_type):
    # calculating air change rate according to the guide lines in SAP 2012 Appendix P
    if (op_type == 0) and (bld_type == 0):
        return float(0.1)
    elif (op_type == 0) and (bld_type == 1):
        return float(0.1)
    elif (op_type == 0) and (bld_type == 2):
        return float(0.2)
    elif (op_type == 0) and (bld_type == 3):
        return float(0.1)
    elif (op_type == 1) and (bld_type == 0):
        return float(0.8)
    elif (op_type == 1) and (bld_type == 1):
        return float(0.5)
    elif (op_type == 1) and (bld_type == 2):
        return float(1)
    elif (op_type == 1) and (bld_type == 3):
        return float(0.6)
    elif (op_type == 2) and (bld_type == 0):
        return float(3)
    elif (op_type == 2) and (bld_type == 1):
        return float(2)
    elif (op_type == 2) and (bld_type == 2):
        return float(4)
    elif (op_type == 2) and (bld_type == 3):
        return float(2.5)
    elif (op_type == 3) and (bld_type == 0):
        return float(6)
    elif (op_type == 3) and (bld_type == 1):
        return float(4)
    elif (op_type == 3) and (bld_type == 2):
        return float(8)
    elif (op_type == 3) and (bld_type == 3):
        return float(5)
    else:
        print("Please enter the building type and opening type according to the guide line")
        
def guid_to_brep(breps):
    # this definition converts the brep list to GUID in order to be used for union process
    list = []
    for brep in breps:
        list.append(doc.Objects.AddBrep(brep))
    return list


def solar_flux_calculator(surface, s_flux, latitude, North, max_roof_angle, month_Solar_declination):
    if not North:
        North = rg.Vector3d.YAxis
    else:
        North = North
    k_values = { 'North' : [26.3,-38.5,14.8,-16.5,27.3 ,-11.9,-1.06, 0.0872, -0.191],
    'South' : [-0.66,-0.106,2.93,3.63,-0.374,-7.4,-2.71,-0.991,4.59],
    'SE' : [-2.95,2.89,1.17,5.67,-3.54,-4.28,-2.72,-0.25,3.07],
    'East' : [1.44,-2.36,1.07,-0.514,1.89,-1.64,-0.542,-0.757,0.604],
    'NE' : [0.165,-3.68,3.0 ,6.38 ,-4.53 ,-0.405 ,-4.38 ,4.89 ,-1.99]}
            
    surface_vector = surface.NormalAt(0.5,0.5)
    delta = math.radians(latitude - month_Solar_declination)
    # calculating surface inclination
    srf_incl = rg.Vector3d.VectorAngle(surface_vector, rg.Vector3d.ZAxis) * 180 / math.pi
            
    # calculating conveting factor
    if srf_incl <= max_roof_angle:
        Rh = 1
        key = 0
    else:
        # calculating the orientation
        vector_angle = rg.Vector3d.VectorAngle(surface_vector, North) * 180 / math.pi
        if vector_angle < 36:
            key = 'North'
        if vector_angle >= 36 and vector_angle < 72:
            key = 'NE'
        if vector_angle >= 72 and vector_angle < 108:
            key = 'East'
        if vector_angle >= 108 and vector_angle < 144:
            key = 'SE'
        if vector_angle >= 144 and vector_angle <= 180:
            key = 'South'
    A = k_values[key][0] * (math.sin(math.radians(srf_incl/2)) ** 3) + k_values[key][1] * ((math.radians(srf_incl/2)) ** 2) + k_values[key][2] * (math.radians(srf_incl/2))
    B = k_values[key][3] * (math.sin(math.radians(srf_incl/2)) ** 3) + k_values[key][4] * ((math.radians(srf_incl/2)) ** 2) + k_values[key][5] * (math.radians(srf_incl/2))
    C = k_values[key][6] * (math.sin(math.radians(srf_incl/2)) ** 3) + k_values[key][7] * ((math.radians(srf_incl/2)) ** 2) + k_values[key][8] * (math.radians(srf_incl/2)) + 1
    Rh = (A * (math.cos(delta)**2)) + (B * (math.cos(delta))) + C
    return (s_flux * Rh, key, srf_incl)


def z_summer(surface, window_height, z_blinds, solar_acces, overhang, time_fraction, orientation, srf_incl, WWR, max_roof_angle):
    if srf_incl < max_roof_angle:
        return 1
    if solar_acces > 80:
         solar_acces_factor = 0.5
    elif solar_acces > 60 and solar_acces <= 80:
        solar_acces_factor = 0.7
    elif solar_acces > 20 and solar_acces <= 60:
        solar_acces_factor = 0.9
    elif solar_acces <= 20:
        solar_acces_factor = 1
    if not overhang and not z_blinds and not time_fraction:
        return solar_acces_factor
    if not window_height:
        window_height = math.sqrt(surface.GetSurfaceSize()[2])
    # defining a dictionary for Shading factors for blinds, curtains or external shutters according to SAP table P3
    blind_dic = { '0' : 0.8 , '1' : 0.9, '2' : 0.85, '3' : 0.6, '4' : 0.88 , '5' : 0.7, '6' : 0.27, '7' : 0.24, '8' : 0.85, '9' : 0.65}
    # Defining the shading factor of the blinds according to the fraction of the daylight hours
    if time_fraction and z_blinds:
        blind_factor = time_fraction * blind_dic ['{}'.format(z_blinds)] + (1 - time_fraction)
    elif not time_fraction and z_blinds:
        blind_factor = blind_dic ['{}'.format(z_blinds)]
    elif not z_blinds:
        blind_factor = 1
    # finding the depth factor of the overhang window 
    d_factor = overhang / window_height
    if d_factor < 0.1:
        depth_factor = 'zr'
    elif d_factor >= 0.1 and d_factor < 0.3:
        depth_factor = 'two'
    elif d_factor >= 0.3 and d_factor < 0.5:
        depth_factor = 'four'
    elif d_factor >= 0.5 and d_factor < 0.7:
        depth_factor = 'six'
    elif d_factor >= 0.7 and d_factor < 0.9:
        depth_factor = 'eight'
    elif d_factor >= 0.9 and d_factor < 1.1:
        depth_factor = 'ten'
    else:
        depth_factor = 'twlv'
    window_width = surface.GetSurfaceSize()[1]
    if overhang/window_width < 2:
        overhang_dic = {'Northzr': 1, 'NEzr': 1, 'Eastzr': 1, 'SEzr': 1, 'Southzr': 1, 
        'Northtwo': 0.94, 'NEtwo': 0.91, 'Easttwo': 0.89, 'SEtwo': 0.84, 'Southtwo': 0.79,
        'Northfour': 0.9, 'NEfour': 0.85, 'Eastfour': 0.79, 'SEfour': 0.72, 'Southfour': 0.64,
        'Northsix': 0.88, 'NEsix': 0.81, 'Eastsix': 0.72, 'SEsix': 0.62, 'Southsix': 0.53,
        'Northeight': 0.86, 'NEeight': 0.79, 'Easteight': 0.66, 'SEeight': 0.55, 'Southeight': 0.5,
        'Northten': 0.85, 'NEten': 0.77, 'Eastten': 0.61, 'SEten': 0.52, 'Southten': 0.49,
        'Northtwlv': 0.84, 'NEtwlv': 0.76, 'Easttwlv': 0.57, 'SEtwlv': 0.5, 'Southtwlv': 0.48}
    else:
        overhang_dic = {'Northzr': 1, 'NEzr': 1, 'Eastzr': 1, 'SEzr': 1, 'Southzr': 1, 
        'Northtwo': 0.92, 'NEtwo': 0.89, 'Easttwo': 0.88, 'SEtwo': 0.83, 'Southtwo': 0.79,
        'Northfour': 0.85, 'NEfour': 0.8, 'Eastfour': 0.76, 'SEfour': 0.67, 'Southfour': 0.55,
        'Northsix': 0.79, 'NEsix': 0.72, 'Eastsix': 0.66, 'SEsix': 0.54, 'Southsix': 0.38,
        'Northeight': 0.73, 'NEeight': 0.65, 'Easteight': 0.58, 'SEeight': 0.43, 'Southeight': 0.32,
        'Northten': 0.69, 'NEten': 0.59, 'Eastten': 0.51, 'SEten': 0.36, 'Southten': 0.3,
        'Northtwlv': 0.66, 'NEtwlv': 0.55, 'Easttwlv': 0.46, 'SEtwlv': 0.31, 'Southtwlv': 0.29}
    return float(blind_factor * (solar_acces_factor+ overhang_dic[orientation + depth_factor] - 1))

def temprature_threshold(TMP, Summer_ratio, ex_t_average):
    # Obtaining the threshold internal temperature which is used to estimate likelihood of high internal temperature.
    if TMP < 285:
        delta_t_mass = 2.0 - 0.007 * TMP
    else:
        delta_t_mass = 0
    temp_thresh = ex_t_average + Summer_ratio + delta_t_mass
    if temp_thresh < 20.5:
        overheating = 'Not significant'
        color = "3,166,90"
    elif temp_thresh >= 20.5 and temp_thresh < 22:
        overheating = 'Slight'
        color = "242,191,39"
    elif temp_thresh >= 22 and temp_thresh < 23.5:
        overheating = 'Medium'
        color = "217,121,4"
    else:
        overheating = 'High'
        color = "191,4,4"
    return overheating, color


# making the union of the breps to exclude the inner partions of the faces
b1 = rs.BooleanUnion(guid_to_brep(masses))

# convering the new unioned brep GUID into brep class
brep_union = rs.coercebrep(b1)

results = []
preview = []

for i, brep in enumerate(masses):
    loss  = 0
    gain = 0
    for face in brep_union.Faces:
        # Calculating the Trimed faces center points in order to find which zone the new face belongs to 
        center_point = rg.AreaMassProperties.Compute(face).Centroid
        # Calculating the heat gain and heat loss of each external surface
        if ghc.PointInBrep(brep, center_point, False) and ((rg.Vector3d.VectorAngle(face.NormalAt(0.5,0.5), rg.Vector3d.ZAxis) * 180 / math.pi) < 170) and ((rg.Vector3d.VectorAngle(face.NormalAt(0.5,0.5), rg.Vector3d.ZAxis) * 180 / math.pi) > max_roof_angle):
            s_flux, key, srf_incl = solar_flux_calculator(face, solar_flux, Latitude, North, max_roof_angle, month_Solar_declination)
            wall_area = rg.AreaMassProperties.Compute(face).Area
            shade_factor = z_summer(face, window_height, curt_blind, solar_access, overhangs, time_fraction, key, srf_incl, WWR, max_roof_angle)
            s_gain = 0.9 * wall_area * WWR * FF * s_flux * shade_factor * g_value
            gain += s_gain
            heat_loss = (wall_area * (1-WWR) * U_wall) + (wall_area * WWR * U_window) + 0.15 * wall_area
            loss += heat_loss
        elif ghc.PointInBrep(brep, center_point, False) and ((rg.Vector3d.VectorAngle(face.NormalAt(0.5,0.5), rg.Vector3d.ZAxis) * 180 / math.pi) <= max_roof_angle) and Roof:
            loss += (rg.AreaMassProperties.Compute(face).Area) * (U_roof+0.15)
        else:
            continue
    # calculating the volume of the space
    vol = rg.VolumeMassProperties.Compute(brep).Volume
    # calculating ventilation heat loss:
    n = air_change_rate(opn_type,bld_type)
    vent_loss = 0.33 * n * vol
    loss += vent_loss
    ratio = gain / loss
    results.append(temprature_threshold(TMP, ratio, ex_t_average)[0])
    # coloring the mesh according to the temprature threshold
    for mesh in rg.Mesh.CreateFromBrep(brep):
        preview.append(ghc.MeshColours(mesh,(temprature_threshold(TMP, ratio, ex_t_average)[1])))