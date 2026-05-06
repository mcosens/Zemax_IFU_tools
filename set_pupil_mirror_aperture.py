'''
set_pupil_mirror_aperture.py
Written by: Maren Cosens
Date: 5/4/26

Description: uses ZOS-API to set the aperture of the IFU pupil mirrors based on the footprints across all MIRMOS bands
'''
##import packages
import os
import pandas as pd
import zos_pyclass
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from System import Enum, Int32, Double
from zemax_functions import get_footprint, local_to_global_coords #would need to make more general and pass api instance to make this work - Todo later

##set-up rcparams to control plot style (later move to a custom file that is read in)
#legend
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.labelspacing'] = 0.2
mpl.rcParams['legend.handletextpad'] = 0.1
mpl.rcParams['legend.columnspacing'] = 1
mpl.rcParams['legend.borderaxespad'] = 0.1
#fontsizes
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['axes.labelsize'] = 12

##main script
# load local variables
zos = zos_pyclass.PythonStandaloneApplication()

ZOSAPI = zos.ZOSAPI
IFU_System = zos.TheSystem
TheApplication = zos.TheApplication
    
# Setup
IFU_System.LoadFile(r'C:\Users\mcosens\Documents\Zemax\MIRMOS\IFU\MIRMOS_full_IFS.zos', False)
#path to save results
res_dir='C:\\Users\\mcosens\\Documents\\Research_docs\\MIRMOS\\IFU\\'

#set to config1 to start
IFU_MCE=IFU_System.MCE
IFU_MCE.SetCurrentConfiguration(1)
#get system paramaters
IFU_LDE=IFU_System.LDE
nfields=IFU_System.SystemData.Fields.NumberOfFields
nwaves=IFU_System.SystemData.Wavelengths.NumberOfWavelengths
nconfigs=IFU_MCE.NumberOfConfigurations
bands=['i', 'Y', 'J', 'H', 'K']
#pupil mirror surface number
pmirr_surf = 41

#get footprints at pupil mirrors 
for c in range(nconfigs):
    get_footprint(ZOSAPI, IFU_System, config=c+1, field=0, wave=0, nrays=10, surface=pmirr_surf, delete_vignetted=False, outpath=res_dir+'footprints\\', outfile=f'footprint_config{c+1}_allfields_allwaves_pupil.txt')
#read in text file and get values for center and radii
pmir_xcen = np.full(nconfigs, fill_value=np.nan)
pmir_ycen = np.full(nconfigs, fill_value=np.nan)
pmir_xrad = np.full(nconfigs, fill_value=np.nan)
pmir_yrad = np.full(nconfigs, fill_value=np.nan)
for c in range(nconfigs):
    with open(f"{res_dir}\\footprints\\footprint_config{c+1}_allfields_allwaves_pupil.txt", 'r', encoding='utf-16') as f:
        lines = f.readlines()
        for line in lines:
            if '=' in(line):
                param, value = line.split('=')
                if param.strip()=='Ray X Center':
                    pmir_xcen[c] = float(value.strip())
                elif param.strip()=='Ray Y Center':
                    pmir_ycen[c] = float(value.strip())
                elif param.strip()=='Ray X Half Width':
                    pmir_xrad[c] = float(value.strip())
                elif param.strip()=='Ray Y Half Width':
                    pmir_yrad[c] = float(value.strip())
## add rectangular aperture to pupil mirror surface
aperture = IFU_LDE.GetSurfaceAt(pmirr_surf).ApertureData.CreateApertureTypeSettings(ZOSAPI.Editors.LDE.SurfaceApertureTypes.RectangularAperture)
aperture._S_RectangularAperture.XHalfWidth = 20 #arbitrary for now, will be set by analysis below
aperture._S_RectangularAperture.YHalfWidth = 20
IFU_LDE.GetSurfaceAt(pmirr_surf).ApertureData.ChangeApertureTypeSettings(aperture)

## add MCE operands to control aperture
new_op_types = [ZOSAPI.Editors.MCE.MultiConfigOperandType.MOFF, ZOSAPI.Editors.MCE.MultiConfigOperandType.APDX, ZOSAPI.Editors.MCE.MultiConfigOperandType.APDY, ZOSAPI.Editors.MCE.MultiConfigOperandType.APMN, ZOSAPI.Editors.MCE.MultiConfigOperandType.APMX] #decenter-x, decenter-y, x half-width, y half-width
## calculate overall centers and x,y extent for a given slice (need to use every 21st entry of arrays)
new_op_params=np.full((len(new_op_types)-1,int(nconfigs/len(bands))), fill_value=np.nan)
for i in range(int(nconfigs/len(bands))):
    indices=i+21*np.arange(len(bands))
    new_op_params[0,i] = np.nanmean(pmir_xcen[indices]) #mean x center
    new_op_params[1,i] = np.nanmean(pmir_ycen[indices]) #mean y center
    new_op_params[2,i] = np.nanmax(pmir_xrad[indices])*1.05 #mirror x half-width (max ray radius - decenter) - add 5% clear aperture
    new_op_params[3,i] = np.nanmax(pmir_yrad[indices])*1.05 #mirror y half-width
n_operands= IFU_MCE.NumberOfOperands #current number of operands
for i in range(len(new_op_types)):
    IFU_MCE.AddOperand()
    current_op=IFU_MCE.GetOperandAt(n_operands+i+1)
    current_op.ChangeType(new_op_types[i])
    if i!=0:
        current_op.Param1=pmirr_surf
        #set values for i-band configs based on all bands and add pickup for rest
        for j in range(1,nconfigs+1):
            if j<22:
                current_op.GetOperandCell(j).Value = str(new_op_params[i-1,j-1])
            else:
                pickup_solve=current_op.GetOperandCell(j).CreateSolveType(ZOSAPI.Editors.SolveType.ConfigPickup)
                pickup_solve._S_ConfigPickup.Configuration = j-int((j-1)/21)*21
                pickup_solve._S_ConfigPickup.Operand = n_operands+i+1
                current_op.GetOperandCell(j).SetSolveData(pickup_solve)
    else:
        current_op.GetOperandCell(1).Value = 'Pupil Mirror Aperture'
        
#check for interference of adjacent mirrors by getting global x,y,z of each mirror

#pupil mirror layout schematic (number = config)
# 11   8    5    2    13    16    19
#   10   7    4    1    14    17     20
#      9   6    3   12     15    18     21

#x-direction interference:
slice_order_x=np.concatenate([np.linspace(11,1,11), np.linspace(12,21,10)])
mirror_space_x=np.full(len(slice_order_x)-3, fill_value=np.nan)
mirror_space_xcen=np.full(len(slice_order_x)-3, fill_value=np.nan)
for i in range(len(mirror_space_x)):
    c = int(slice_order_x[i])
    c_adj = int(slice_order_x[i+3])
    #get local coordinates,aperture size and rotation matrix for each and convert to global coordinates
    #current congfig
    IFU_MCE.SetCurrentConfiguration(c)
    rot_matrix_c = IFU_LDE.GetGlobalMatrix(pmirr_surf)
    #use aperture coordinates for local coords
    lefttop_edge_loc = [new_op_params[0,c-1]-new_op_params[2,c-1], new_op_params[1,c-1]+new_op_params[3,c-1],0] #[x,y,z]
    leftbot_edge_loc = [new_op_params[0,c-1]-new_op_params[2,c-1], new_op_params[1,c-1]-new_op_params[3,c-1],0] #[x,y,z]
    leftcen_loc = [new_op_params[0,c-1]-new_op_params[2,c-1], new_op_params[1,c-1],0] #[x,y,z]
    lefttop_edge = coordinate_transform(rot_matrix_c[1:], lefttop_edge_loc)
    leftbot_edge = coordinate_transform(rot_matrix_c[1:], leftbot_edge_loc)
    left_cen =  coordinate_transform(rot_matrix_c[1:], leftcen_loc)
    left_extreme=np.min((leftbot_edge[0], lefttop_edge[0]))
    #adjacent configuration
    IFU_MCE.SetCurrentConfiguration(c_adj)
    rot_matrix_c_adj = IFU_LDE.GetGlobalMatrix(pmirr_surf)
    righttop_edge_loc = [new_op_params[0,c_adj-1]+new_op_params[2,c_adj-1], new_op_params[1,c_adj-1]+new_op_params[3,c_adj-1],0] #[x,y,z]
    rightbot_edge_loc = [new_op_params[0,c_adj-1]+new_op_params[2,c_adj-1], new_op_params[1,c_adj-1]-new_op_params[3,c_adj-1],0] #[x,y,z]
    rightcen_loc = [new_op_params[0,c_adj-1]+new_op_params[2,c_adj-1], new_op_params[1,c_adj-1], 0] #[x,y,z]
    righttop_edge = coordinate_transform(rot_matrix_c_adj[1:], righttop_edge_loc)
    rightbot_edge = coordinate_transform(rot_matrix_c_adj[1:], rightbot_edge_loc)
    right_cen = coordinate_transform(rot_matrix_c_adj[1:], rightcen_loc)
    right_extreme=np.max((rightbot_edge[0], righttop_edge[0]))
    #calculate/save space between mirrors (negative values indicate overlap)
    mirror_space_x[i]=left_extreme-right_extreme
    mirror_space_xcen[i] = left_cen[0] - right_cen[0]

#y-direction interference:
slice_order_y = [11, 8, 5, 2, 13, 16, 19, 10, 7, 4, 1, 14, 17, 20]
mirror_space_y = np.full(len(slice_order_y), fill_value=np.nan)
cen_x_offset = np.full(len(slice_order_y), fill_value=np.nan)
for i in range(len(mirror_space_y)):
    c = int(slice_order_y[i])
    if c==1:
        c_adj=12
    elif c < 13:
        c_adj = c-1
    else:
        c_adj = c+1
    #get local coordinates,aperture size and rotation matrix for each and convert to global coordinates
    #current congfig
    IFU_MCE.SetCurrentConfiguration(c)
    rot_matrix_c = IFU_LDE.GetGlobalMatrix(pmirr_surf)
    #use aperture coordinates for local coords
    topcen_loc = [new_op_params[0,c-1], new_op_params[1,c-1]+new_op_params[3,c-1],0] #[x,y,z]
    top_cen =  coordinate_transform(rot_matrix_c[1:], topcen_loc)
    #adjacent configuration
    IFU_MCE.SetCurrentConfiguration(c_adj)
    rot_matrix_c_adj = IFU_LDE.GetGlobalMatrix(pmirr_surf)
    botcen_loc = [new_op_params[0,c_adj-1], new_op_params[1,c_adj-1]-new_op_params[3,c_adj-1], 0] #[x,y,z]
    bot_cen = coordinate_transform(rot_matrix_c_adj[1:], botcen_loc)
    #calculate/save space between mirrors (negative values indicate overlap)
    mirror_space_y[i] = top_cen[1] - bot_cen[1]
    cen_x_offset[i] = np.abs(top_cen[0] - bot_cen[0])
mirror_center_offset_r = np.sqrt(mirror_space_y**2 + cen_x_offset**2)

#save updated model
IFU_System.SaveAs("C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\MIRMOS_full_IFS_apertures.zos")  #save as seperate file to test
# close server instance of OpticStudio
del zos
zos = None