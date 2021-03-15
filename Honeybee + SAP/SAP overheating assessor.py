"""Use this plugin to take any list of Honeybee rooms to calculate the risk of overheating. The component uses the properties of the honeybee rooms to calculate room's propeties. Singel_sided and double_sided ventilation conditions are calculated according the window's operability.
    Inputs:
        north_: A vector to represent the North direction. If not provided the default Y axis will be selected as North
        _latitude: Latitude of the locatation of the building
        _hb_objs: A list of honeybee rooms to perform the overheating study on the building masses
        _bld_type: Building type according to the SAP Table p1. The values: {"single story dwelling with possibility of cross ventilation" = 0 , "single story dwelling cross ventilation not possible" = 1, "Dwelling of two or more storeys windows open upstairs and downstairs Cross ventilation possible" = 2, "Dwelling of two or more storeys windows open upstairs and downstairs Cross ventilation not possible" =3} 
        _solar_acces_: Percentage of Sky blocked by obstecles. {"Heavy" > 80, "More than avarage" >60 & <80, "average or Unknown" >20 & <60, "very little" <20} 
        _solar_flux: solar flux on the horizontal surface during the summer period
        z_blinds_: Selcet the type of blinds, curtains or external shutters to calculate the shading factor according to SAP table P3: {"Net curtain (covering whole window)" = 0, " Net curtain (covering half window)" = 1, " Dark-coloured curtain or roller blind (note 1)" = 2, " Light-coloured curtain or roller blind (note 1)" = 3, " Dark-coloured venetian blind (note 2)" = 4," Light-coloured venetian blind (note 2)” = 5, “Dark-coloured external shutter, window closed (notes 1, 3)" = 6, “White external shutter, window closed (notes 1, 3)” = 7, “Dark-coloured external shutter, window fully open (notes 1, 3)” = 8, “White external shutter, window fully open (notes 1, 3)” = 9}
        _op_type: Ventilation Opening type according to the SAP Table p1.he values: {"Trickle vents only" = 0 , "Windows slightly open (50 mm)" = 1, "Windows open half the time" = 2,"Windows fully open" = 3}
        _TMP_: Thermal mass parameter, Defult set to 220
        FF_: The frame factor for windows and doors, Defult "0.7" (fraction of opening that is glazed). Frame factor according to material: {"wood" = 0.7, "metal" = 0.8, PVC-U = 0.7} 
        _ex_t_average: The mean external temperature for the selected summer month and climate region.
        time_fraction_: fraction of the daylight hours that shading blinds are in use
        month_: month in which you want to calculate overheating risk for, Defult set for July. {June : "0", July :"1", August : "2"} 
    Outputs:
        results: Likelihood of high internal temperature during hot weather.
        preview: Colored mesh preview
"""
#Copyright (c) 2021, Hamidreza Shahriari 
__author__ = "Hamidreza"
import Rhino.Geometry as rg
from ladybug_rhino.togeometry import to_vector2d, to_point2d, to_point3d
from ladybug_geometry.geometry2d.pointvector import Vector2D, Point2D
import math
from ladybug_geometry.geometry3d.pointvector import Vector3D
from honeybee.facetype import Wall, RoofCeiling, Floor
from honeybee.boundarycondition import Outdoors
import ladybug.color as lc
from ladybug_rhino.fromgeometry import from_face3ds_to_colored_mesh, from_face3d_to_wireframe

if north_:
    north_ = to_vector2d(north_)
else:
    north_ = to_vector2d(rg.Vector3d(0,1,0))



def solar_gain(face, north, solar_flux , z_blinds, solar_acces, time_fraction, latitude, month , FF):
    s_gain = 0
    s_loss = 0
    for ap in face.apertures:
        winArea = ap.area
        u_val = ap.properties.energy.construction.u_value
        g_val = ap.properties.energy.construction.solar_transmittance
        normal = ap.normal
        dir = ap.cardinal_direction(north)
        # Here we record the directions that operable windows are facing.
        # (Later we use this dictonary to check if the zone is single sided ventelated or double sided)
        if ap.is_operable:
            if dir not in win_dir.keys():
                win_dir[dir] = 1
        s_flux, orientation = solar_flux_calculator(ap, north,solar_flux,latitude, month = '1')
        shade_factor = z_summer(ap, z_blinds, time_fraction, orientation, solar_acces)
        s_loss = winArea * u_val
        if not FF:
            FF = 0.7
        s_gain += 0.9 * winArea * FF * s_flux * shade_factor * g_val
    return (s_gain, s_loss)

def solar_flux_calculator(surface, north,solar_flux,latitude, month ):
    if not month:
        month = "1"
    month_dic = {'0' : 23.1, '1' : 21.2, '2' : 13.7}
    month_Solar_declination = month_dic[month]
    k_values = { 'North' : [26.3,-38.5,14.8,-16.5,27.3 ,-11.9,-1.06, 0.0872, -0.191],
    'South' : [-0.66,-0.106,2.93,3.63,-0.374,-7.4,-2.71,-0.991,4.59],
    'SE' : [-2.95,2.89,1.17,5.67,-3.54,-4.28,-2.72,-0.25,3.07],
    'East' : [1.44,-2.36,1.07,-0.514,1.89,-1.64,-0.542,-0.757,0.604],
    'NE' : [0.165,-3.68,3.0 ,6.38 ,-4.53 ,-0.405 ,-4.38 ,4.89 ,-1.99]}
    dir = face.cardinal_direction(north)
    normal = face.normal
    delta = math.radians(latitude - month_Solar_declination)
    # calculating surface inclination
    srf_incl = normal.angle(Vector3D(0,0,1)) * 180 / math.pi
    if dir == 'NorthEast' or dir == 'NorthWest':
        key = 'NE'
    elif dir == 'SouthWest' or dir == 'SouthEast':
        key = 'SE'
    elif dir == 'East' or dir == 'West':
        key = 'East'
    else:
        key = dir
    A = k_values[key][0] * (math.sin(math.radians(srf_incl/2)) ** 3) + k_values[key][1] * ((math.radians(srf_incl/2)) ** 2) + k_values[key][2] * (math.radians(srf_incl/2))
    B = k_values[key][3] * (math.sin(math.radians(srf_incl/2)) ** 3) + k_values[key][4] * ((math.radians(srf_incl/2)) ** 2) + k_values[key][5] * (math.radians(srf_incl/2))
    C = k_values[key][6] * (math.sin(math.radians(srf_incl/2)) ** 3) + k_values[key][7] * ((math.radians(srf_incl/2)) ** 2) + k_values[key][8] * (math.radians(srf_incl/2)) + 1
    Rh = (A * (math.cos(delta)**2)) + (B * (math.cos(delta))) + C
    return (solar_flux * Rh, key)

def z_summer(surface, z_blinds , time_fraction, orientation , solar_acces = 1):

    if solar_acces > 80:
         solar_acces_factor = 0.5
    elif solar_acces > 60 and solar_acces <= 80:
        solar_acces_factor = 0.7
    elif solar_acces > 20 and solar_acces <= 60:
        solar_acces_factor = 0.9
    elif solar_acces <= 20:
        solar_acces_factor = 1
    # Calculating the overhang depth
    overhang_result = 0
    try:
        shade = surface.shades[0]
        shdBondGmtry = shade.geometry.boundary
        overhang =  shdBondGmtry[0].distance_to_point(shdBondGmtry[3])
        window_height =  surface.geometry.boundary[0].distance_to_point(surface.geometry.boundary[3])
        d_factor = overhang / window_height
        # finding the depth factor of the overhang window 
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
        window_width = surface.geometry.boundary[0].distance_to_point(surface.geometry.boundary[1])
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
        overhang_result += overhang_dic[orientation + depth_factor]
    except:
        if not z_blinds and not time_fraction:
            return solar_acces_factor
    if overhang_result == 0:
        overhang_result = 1
    # defining a dictionary for Shading factors for blinds, curtains or external shutters according to SAP table P3
    blind_dic = { '0' : 0.8 , '1' : 0.9, '2' : 0.85, '3' : 0.6, '4' : 0.88 , '5' : 0.7, '6' : 0.27, '7' : 0.24, '8' : 0.85, '9' : 0.65}
    # Defining the shading factor of the blinds according to the fraction of the daylight hours
    if time_fraction and z_blinds:
        blind_factor = time_fraction * blind_dic ['{}'.format(z_blinds)] + (1 - time_fraction)
    elif not time_fraction and z_blinds:
        blind_factor = blind_dic ['{}'.format(z_blinds)]
    elif not z_blinds:
        blind_factor = 1
    return float(blind_factor * (solar_acces_factor+ overhang_result - 1))

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
    elif (bld_type == -1):
        return float(0.01)
    else:
        print("Please enter the building type and opening type according to the guide line")

def temprature_threshold(TMP, Summer_ratio, ex_t_average):
    if not TMP:
        TMP = 220
    # Obtaining the threshold internal temperature which is used to estimate likelihood of high internal temperature.
    if TMP < 285:
        delta_t_mass = 2.0 - 0.007 * TMP
    else:
        delta_t_mass = 0
    temp_thresh = ex_t_average + Summer_ratio + delta_t_mass
    if temp_thresh < 20.5:
        overheating = 'Not significant'
        color = lc.Color(3,166,90)
    elif temp_thresh >= 20.5 and temp_thresh < 22:
        overheating = 'Slight'
        color = lc.Color(242,191,39)
    elif temp_thresh >= 22 and temp_thresh < 23.5:
        overheating = 'Medium'
        color = lc.Color(217,121,4)
    else:
        overheating = 'High'
        color =  lc.Color(191,4,4)
    return overheating, color

floors = []
ratio = []
wire_view = []

for obj in _hb_objs:
    loss = 0
    gain = 0
    #dictionary to check the direction windows facing for checking single_sided or double_sided vantilation
    win_dir = {}
    for face in obj.faces:
        wire_view.append(from_face3d_to_wireframe(face.geometry))
        bc = face.boundary_condition
        type = face.type
        if isinstance(type, Wall):
            if isinstance(bc, Outdoors):
                s_flux , srf_incl = solar_flux_calculator(face, north_ ,_solar_flux ,_latitude, month_)
                loss += (face.area - face.aperture_area) * face.properties.energy.construction.u_value
                loss += 0.15 * face.area
                win_gain, win_loss = solar_gain(face, north_,_solar_flux, z_blinds_, _solar_acces_, time_fraction_, _latitude, month_, FF_ )
                loss += win_loss
                gain += win_gain
        # Calculating heat loss from the roof
        elif isinstance(type, RoofCeiling):
            if isinstance(bc, Outdoors):
                loss += face.area * face.properties.energy.construction.u_value
                loss += 0.15 * face.area
        elif isinstance(type, Floor):
            floors.append(face)
    if len(win_dir.keys()) > 1 and _bld_type == 0:
        blding_type = 0
    elif len(win_dir.keys()) == 1 and _bld_type == 0:
        blding_type = 1
    elif len(win_dir.keys()) > 1 and _bld_type == 1:
        blding_type = 2
    elif len(win_dir.keys()) == 1 and _bld_type == 1:
        blding_type = 3
    elif len(win_dir.keys()) == 0:
        ratio.append(gain/loss)
        continue
    else:
        print("Please enter a valid number according to the description")
    # calculating heat loss due ventilation
    n = air_change_rate(_op_type,blding_type)
    ven_loss = 0.33 * n * obj.volume
    loss += ven_loss
    ratio.append(gain/loss)

result = []
preview_result = []
a = []

for i in ratio:
    a.append(temprature_threshold(_TMP_, i, _ex_t_average)[1])


for i in range(len(ratio)):
    preview_result.append(from_face3ds_to_colored_mesh([floors[i].geometry], temprature_threshold(_TMP_, ratio[i], _ex_t_average)[1]))
    result.append(temprature_threshold(_TMP_, ratio[i], _ex_t_average)[0])
