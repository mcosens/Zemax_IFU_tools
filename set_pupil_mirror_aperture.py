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
#from make_IFU_plots import get_footprint #would need to make more general and pass api instance to make this work - Todo later

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

    
##functions to be used
def get_footprint(system, config=1, field=0, wave=0, nrays=10, surface=0, delete_vignetted=True, outpath=".\\", outfile="footprint_res.txt"):
    """
    Computes footprint diagram data for given configuration, field #, and wavelength # using ZOS-API analysis.

    Parameters:
    - system: The ZOS system object
    - config (default=1): Configuration number
    - field (default=0): Field number (0=all fields)
    - wave (default=0): Wavelength number (0=all wavelengths)
    - nrays (default=10): Number of rays for the analysis
    - surface (default=0): Surface number to evaluate footprint at (0=image surface)
    - delete_vignetted (default=True): Whether to delete vignetted rays
    - outpath (default=".\\"): filepath to save intermediate and final files
    - outfile (default="footprint_res.txt"): filename to save footprint data (e.g. 'footprint_res.txt')

    Returns:
    - saves footprint diagram data to text file at outfile
    """
    system.MCE.SetCurrentConfiguration(config) #set to first config of K band
    footprint = system.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.FootprintSettings) #Footprint Diagram not yet built in
    #modify settings
    foot_settings = footprint.GetSettings()
    foot_cfgFile = f"{outpath}footprint_settings.cfg"
    foot_settings.SaveTo(foot_cfgFile)
    foot_settings.ModifySettings(foot_cfgFile, "FOO_RAYDENSITY", "10")
    if surface==0:
        foot_settings.ModifySettings(foot_cfgFile, "FOO_SURFACE", str(system.LDE.NumberOfSurfaces))
    else:
        foot_settings.ModifySettings(foot_cfgFile, "FOO_SURFACE", str(surface))
    foot_settings.ModifySettings(foot_cfgFile, "FOO_FIELD", str(field))
    foot_settings.ModifySettings(foot_cfgFile, "FOO_WAVELENGTH", str(wave))
    if delete_vignetted:
        foot_settings.ModifySettings(foot_cfgFile, "FOO_DELETEVIGNETTED", "1") #1=yes, 0=no
    else:
        foot_settings.ModifySettings(foot_cfgFile, "FOO_DELETEVIGNETTED", "0") #1=yes, 0=no
    foot_settings.LoadFrom(foot_cfgFile)
    footprint.ApplyAndWaitForCompletion()
    #read results and save to text file (later can read in and make plots)
    foot_results = footprint.GetResults()
    foot_results.GetTextFile(f"{outpath}\\{outfile}")
    os.remove(foot_cfgFile) #clean up intermediate file
    footprint.Close() #close analysis to avoid errors with limits on total number of analyses

def coordinate_transform(global_matrix, local_coords):
    '''
    Transforms the specified local coordinates to the global coordinate system based on the matrix given
    Parameters:
    global_matrix (len=12 array): [R11, R12, R13, R21, R22, R23, R31, R32, R33, X0, Y0, Z0]
    local_coords (len=2 array): [x,y,z] coordinates in local frame to be transformed with the global_matrix
    Returns (float):
    global_x, global_y, global_z
    '''
    [R11, R12, R13, R21, R22, R23, R31, R32, R33, X0, Y0, Z0] = global_matrix
    local_x, local_y, local_z = local_coords
    global_x = R11*local_x + R12*local_y + R13*local_z + X0 
    global_y = R21*local_x + R22*local_y + R23*local_z + Y0
    global_z = R31*local_x + R32*local_y + R33*local_z + Z0
    return global_x, global_y, global_z


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
    get_footprint(IFU_System, config=c+1, field=0, wave=0, nrays=10, surface=pmirr_surf, delete_vignetted=False, outpath=res_dir+'footprints\\', outfile=f'footprint_config{c+1}_allfields_allwaves_pupil.txt')
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