'''
tolerance_IFU.py
Written by: Maren Cosens
Date: 5/12/26

Description: uses ZOS-API to perform a tolerance analysis on each slice of the MIRMOS IFU (K band only to start)
1. Use saved tolerance data editor parameters for IFU surfaces / starting points
2. Perform inverse tolerances one configuration at a time to get tightest tolerances for all (likely easiest to apply accross the board)
3. Run MC sensitivity analysis on each configuration
'''
import os
import pandas as pd
import zos_pyclass
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.labelspacing'] = 0.2
mpl.rcParams['legend.handletextpad'] = 0.1
mpl.rcParams['legend.columnspacing'] = 1
mpl.rcParams['legend.borderaxespad'] = 0.1
#fontsizes
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['axes.labelsize'] = 12

# Setup
zos = zos_pyclass.PythonStandaloneApplication()
ZOSAPI = zos.ZOSAPI
IFU_System = zos.TheSystem
IFU_System.LoadFile(r'C:\Users\mcosens\Documents\Zemax\MIRMOS\IFU\virtualIFU_Kband_4-2-26_tol.zmx', False) #includes coordinate break for compensator on doublet position (tip/tilt and decenter)

tol_file = r'C:\Users\mcosens\Documents\Zemax\MIRMOS\IFU\Tolerance\IFU_precision_wcomp.TOL'
#path to save results
res_dir='C:\\Users\\mcosens\\Documents\\Research_docs\\MIRMOS\\IFU\\'

#set to config1 to start
IFU_MCE=IFU_System.MCE
IFU_MCE.SetCurrentConfiguration(1)
#get system paramaters
nfields=IFU_System.SystemData.Fields.NumberOfFields
nwaves=IFU_System.SystemData.Wavelengths.NumberOfWavelengths
nconfigs=IFU_MCE.NumberOfConfigurations

##first get single configuration files set up
for i in range(1,nconfigs+1):
    single_IFU=IFU_System.CopySystem() #make copy of file
    single_IFU.MCE.SetCurrentConfiguration(i)
    single_IFU.MCE.MakeSingleConfiguration() #remove all but current configuration
    single_IFU.Tools.RemoveAllVariables()#remove all variables 
    #remove solves using design lockdown tool (from forum answer: https://community.zemax.com/zpl-13/delete-solves-in-zemax-using-zpl-3325)
    design_lkd = single_IFU.Tools.OpenDesignLockdown()
    design_lkd.ExcludePickups = False
    design_lkd.UsePrecisionRounding = False
    design_lkd.FixModelGlass = False
    design_lkd.ConvertSDToMaxApertures = False
    design_lkd.RunAndWaitForCompletion()
    #load tolerance parameters
    single_IFU.TDE.LoadToleranceFile(tol_file)
    '''
    #remove grating to eliminate spot radius growth due to dispersion
    #ignore surfaces 65-68, remove tilt_x on 67, remove decenter_y on 68
    single_IFU.LDE.RemoveSurfacesAt(65,4) #removes 4 surfaces starting at 63
    s69=single_IFU.LDE.GetSurfaceAt(65) #was surface 69
    s69.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par3).DoubleValue = 0.0
    s70=single_IFU.LDE.GetSurfaceAt(66)
    s70.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par2).DoubleValue = 0.0
    '''
    #instead of above, simply change to single wavelength since there is lateral color induced anyway which we don't care about in this tolerancing step
    single_IFU.SystemData.Wavelengths.RemoveWavelength(3)
    single_IFU.SystemData.Wavelengths.RemoveWavelength(1)
    single_IFU.SaveAs(f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config{i}.zos") 


#close API instance, need a new one for each configuration to evaluate
del zos
zos = None
for i in range(1, nconfigs+1):
    #setup
    zos = zos_pyclass.PythonStandaloneApplication()
    ZOSAPI = zos.ZOSAPI
    single_IFU = zos.TheSystem
    single_IFU.LoadFile(f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config{i}.zos", False)
    #inverse tolerances, save values to get overall minima (increment=0.005)
    single_tol=single_IFU.Tools.OpenTolerancing()
    single_tol.SetupMode = ZOSAPI.Tools.Tolerancing.SetupModes.InverseIncrement
    single_tol.SetupChange = ZOSAPI.Tools.Tolerancing.SetupChanges.LinearDifference
    single_tol.MaximumCriteria = 0.0005 #was 1um, now .5 increase per tolerance operand as a starting point
    # Select Criterion and related settings
    single_tol.Criterion = ZOSAPI.Tools.Tolerancing.Criterions.RMSSpotRadius
    single_tol.CriterionSampling = 3
    single_tol.CriterionField = ZOSAPI.Tools.Tolerancing.CriterionFields.UserDefined
    single_tol.CriterionComp = ZOSAPI.Tools.Tolerancing.CriterionComps.OptimizeAll_DLS #index=3
    # Select number of MC runs and files to save
    single_tol.NumberOfRuns = 20
    single_tol.NumberToSave = 0
    #file to save to
    single_tol.TolDataFile= (f"Kband_config{i}.ZTD") #cannot include path or this will fail; must always save at same location as lens file
    single_tol.OutputFile = f"Kband_config{i}.txt" #try saving as plaintext file; can read from this but not that well formatted, may be easier to use API functions
    single_tol.SaveTolDataFile = True
    # Run the Tolerancing analysis
    single_tol.RunAndWaitForCompletion()
    single_tol.Close()
    #save file
    single_IFU.SaveAs(f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config{i}_sens.zos") 
    del zos
    zos=None

##inspect results of each file 
#generate histogram of impact of each tolerance parameter and illustrate what the system is most sensitive to
#use to determine which tolerances can and should be tightened
col_list = ['Parameter']
col_list.extend(['Change_min_'+str(i) for i in range(1, nconfigs+1)]) #effect on RMS spot radius
col_list.extend(['Change_max_'+str(i) for i in range(1, nconfigs+1)])
col_list.extend(['Tol_min_'+str(i) for i in range(1, nconfigs+1)]) #saved tolerance for each parameter
col_list.extend(['Tol_max_'+str(i) for i in range(1, nconfigs+1)])
df_sens=pd.DataFrame(columns=col_list, index=np.arange(1, nconfigs+1)) #dataframe to hold sensitivity results
for i in range(1, nconfigs+1):
    #setup
    zos = zos_pyclass.PythonStandaloneApplication()
    ZOSAPI = zos.ZOSAPI
    single_IFU = zos.TheSystem
    single_IFU.LoadFile(f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config{i}_sens.zos", False)
    #load tolerance data file to get resulting change in spot radius for each operand
    TolReader=single_IFU.Tools.OpenToleranceDataViewer()
    TolReader.FileName= f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config{i}.ZTD"
    TolReader.RunAndWaitForCompletion()
    sdata = TolReader.SensitivityData
    n_operands = sdata.NumberOfResultOperands
    for id in range(n_operands):
        op = sdata.GetOperand(id)
        op_name=str(op.OperandType) #get name and surface
        for j in range(op.NumberOfParameters):
            if str(op.GetParameter(j).Name) == 'Surface' or str(op.GetParameter(j).Name) == 'Surface1' or str(op.GetParameter(j).Name) == 'Surface2':
                op_name+=f"_s{op.GetParameter(j).IntValue}"
            elif str(op.GetParameter(j).Name) == 'Adjust':
                op_name+=f"_a{op.GetParameter(j).IntValue}"
        change = op.GetEffectOnCriterion(0)
        change_max = change.EstimatedChangeMaximum
        change_min = change.EstimatedChangeMinimum
        if i == 1:
            df_sens.loc[id, 'Parameter'] = op_name #only need to set this once since it is the same for all configurations
        df_sens.loc[id, f'Change_min_{i}'] = change_min
        df_sens.loc[id, f'Change_max_{i}'] = change_max
    TolReader.Close()
    #open tolerance data editor to get resulting tolerace for each operand
    single_TDE=single_IFU.TDE
    n_operands=single_TDE.NumberOfOperands
    skip_par = ['TWAV', 'COMP', 'CPAR'] #parameter names to skip (test wavelength, compensators)
    for id in range(1,n_operands+1):
        op=single_TDE.GetOperandAt(id)
        op_name=str(op.Type) #get name
        if op_name in(skip_par):
            continue
        par1, par2 = op.Param1, op.Param2
        tol_min, tol_max = op.Min, op.Max
        if op.Param1Cell.Header == 'Surf' or op.Param1Cell.Header == 'Surf1':
            op_name+=f"_s{par1}"
        if op.Param2Cell.Header == 'Adj':
            op_name+=f"_a{par2}"
        elif op.Param2Cell.Header == 'Surf2':
            op_name+=f"_s{par2}"
        df_sens.loc[df_sens['Parameter']==op_name, f'Tol_min_{i}'] = tol_min
        df_sens.loc[df_sens['Parameter']==op_name, f'Tol_max_{i}'] = tol_max
    del zos
    zos=None
    
##summarize results accross all configurations
df_summary=pd.DataFrame(columns=['Parameter', 'mean_change', 'max_change', 'min_change', 'low_tol', 'high_tol'], index=np.arange(df_sens.shape[0]))
change_index = int((len(df_sens.columns)-1)/2) #number of change columns for each parameter
for i in range(df_sens.shape[0]):
    mean_change = df_sens.loc[i].values[1:change_index+1].mean()
    max_change = df_sens.loc[i].values[1:change_index+1].max()
    min_change = df_sens.loc[i].values[1:change_index+1].min()
    l_tol = df_sens.loc[i].values[change_index+1:change_index+1+nconfigs].max() #these will be negative to need max
    h_tol = df_sens.loc[i].values[change_index+1+nconfigs:].min() #these will be positive so need min
    df_summary.loc[i] = [df_sens.loc[i, 'Parameter'], mean_change, max_change, min_change, l_tol, h_tol]
#sort by mean change
df_summary.sort_values(by='mean_change', ascending=False, inplace=True)
#print top 10 most sensitive parameters
print(df_summary.head(10))
#most sensitive to tilt and decenter

##plot histogram of sensitivity results / determine worst offenders across all configurations
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.gist_rainbow(np.linspace(0,1,df_sens.shape[0])))
fig=plt.figure(figsize=(10,6))
for i in range(df_sens.shape[0]):
    plt.hist(df_sens.loc[i].values[1:change_index+1]*1000, alpha=0.3, label=str(df_sens.loc[i, 'Parameter']))
plt.xlabel('RMS Spot Radius Change [um]')
plt.title('Sensitivity Analysis')
plt.legend(ncols=3, fontsize='small')
plt.tight_layout()
plt.savefig(f"{res_dir}sensitivity_histogram.png", bbox_inches='tight')
plt.show()

##update tolerance file based on tightest tolerances across all configurations, run longer MC analysis to determine final spot radius
#open config 1 file, update tolerances, save tol file
final_tol_file = r'C:\Users\mcosens\Documents\Zemax\MIRMOS\IFU\Tolerance\IFU_final_tol.TOL'
zos = zos_pyclass.PythonStandaloneApplication()
ZOSAPI = zos.ZOSAPI
single_IFU = zos.TheSystem
single_IFU.LoadFile(f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config1.zos", False)
single_IFU.TDE.LoadToleranceFile(tol_file)
skip_par = ['TWAV', 'COMP', 'CPAR', 'TOFF']
#set new tolerance values
for i in range(1,single_IFU.TDE.NumberOfOperands+1):
    op = single_IFU.TDE.GetOperandAt(i)
    op_name=str(op.Type)
    if op_name in(skip_par):
        continue
    par1, par2 = op.Param1, op.Param2
    if op.Param1Cell.Header == 'Surf' or op.Param1Cell.Header == 'Surf1':
        par_name=op_name+f"_s{par1}"
    if op.Param2Cell.Header == 'Adj':
        par_name+=f"_a{par2}"
    elif op.Param2Cell.Header == 'Surf2':
        par_name+=f"_s{par2}"
    op.Min = df_summary.loc[df_summary['Parameter']==par_name, 'low_tol'].values[0]
    op.Max = df_summary.loc[df_summary['Parameter']==par_name, 'high_tol'].values[0]
#save updated tolerance file
single_IFU.TDE.SaveToleranceFile(final_tol_file)
del zos
zos=None


##get Monte Carlo data (combine with running MC analysis)
df_mc=pd.DataFrame(columns=['RMS_nominal', 'RMS_mean', 'RMS_std', 'RMS_98', 'RMS_90', 'RMS_80', 'RMS_50'], index=np.arange(1, nconfigs+1)) 
CPAR_params={1:'DECX', 2:'DECY', 3:'TETX', 4:'TETY'}
slice_cen=np.zeros(nconfigs)

def run_MC_tol(config, iterations):
    '''Function to run Monte Carlo tolerancing
    Inputs:
    1. config (int): configuration number
    2. iterations (int): number of MC iterations to run
    Returns:
    Run MC and save results to file
    '''
    zos = zos_pyclass.PythonStandaloneApplication()
    ZOSAPI = zos.ZOSAPI
    single_IFU = zos.TheSystem
    single_IFU.LoadFile(f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config{i}.zos", False)
    slice_cen[i-1] = single_IFU.SystemData.Fields.GetField(1).Y
    #load tolerance data file
    single_IFU.TDE.LoadToleranceFile(final_tol_file)
    single_tol=single_IFU.Tools.OpenTolerancing()
    single_tol.SetupMode = ZOSAPI.Tools.Tolerancing.SetupModes.Sensitivity
    single_tol.SetupChange = ZOSAPI.Tools.Tolerancing.SetupChanges.RSSDifference
    # Select Criterion and related settings
    single_tol.Criterion = ZOSAPI.Tools.Tolerancing.Criterions.RMSSpotRadius
    single_tol.CriterionSampling = 3
    single_tol.CriterionField = ZOSAPI.Tools.Tolerancing.CriterionFields.UserDefined
    single_tol.CriterionComp = ZOSAPI.Tools.Tolerancing.CriterionComps.OptimizeAll_DLS #index=3
    # Select number of MC runs and files to save
    single_tol.NumberOfRuns = iterations
    single_tol.NumberToSave = 0
    #file to save to
    single_tol.TolDataFile= (f"Kband_config{i}_MC_{iterations}.ZTD") #cannot include path or this will fail; must always save at same location as lens file
    single_tol.OutputFile = f"Kband_config{i}_MC_{iterations}.txt"
    single_tol.SaveTolDataFile = True
    # Run the Tolerancing analysis
    single_tol.RunAndWaitForCompletion()
    single_tol.Close()
    del zos
    zos=None

#function to get saved MC results without re-running
def get_MC_results(config, df_mc, iterations):
    '''Function to retrieve results of previous Monte Carlo run
    Inputs:
    1. config (int): configuration number
    2. df_mc: DataFrame to store results
    3. iterations (int): number of iterations for run of interest
    Returns:
    Fills in df_mc for row config
    '''
    i=config
    df_mc["iterations"] = iterations
    zos = zos_pyclass.PythonStandaloneApplication()
    ZOSAPI = zos.ZOSAPI
    single_IFU = zos.TheSystem
    single_IFU.LoadFile(f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config{i}.zos", False)
    slice_cen[i-1] = single_IFU.SystemData.Fields.GetField(1).Y
    #load tolerance data file
    #get MC values
    TolReader=single_IFU.Tools.OpenToleranceDataViewer()
    TolReader.FileName= f"C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\Tolerance\\Kband_config{i}_MC_{iterations}.ZTD"
    TolReader.RunAndWaitForCompletion()
    ##Monte Carlo data
    mcData = np.asarray(TolReader.MonteCarloData.Values.Data)
    num_rows, num_cols = mcData.shape
    # Get column name & summary (`GetMetadata` is 0-indexed)
    #need to change to saving desired parts of metadata to list to match with mcData
    for idx in range(num_cols):
        col = TolReader.MonteCarloData.GetMetadata(idx)
        comp_name=None
        if not col.IsOperand:
            if str(col.Name) == 'RmsSpotRadius':
                df_mc.loc[i, 'RMS_nominal'] = col.GetParameter(2).DoubleValue
                df_mc.loc[i, 'RMS_mean'] = np.nanmean(mcData[:,idx]) #same as col.SummaryStatistics.Mean
                df_mc.loc[i, 'RMS_std'] = np.nanstd(mcData[:,idx]) #same as col.SummaryStatistics.PopulationStandardDeviation
                df_mc.loc[i, 'RMS_98'], df_mc.loc[i, 'RMS_90'], df_mc.loc[i, 'RMS_80'], df_mc.loc[i, 'RMS_50'] = np.percentile(mcData[:,idx], q=[98, 90, 80, 50]) #can't get from SummaryStatistics so need full array
            else: #compensators
                n_par=col.NumberOfParameters
                par_list=[]
                for n in range(n_par):
                    par=col.GetParameter(n)
                    par_list.append(str(par.Name))
                if 'Code' in par_list: #COMP
                    comp_name = f"THIC_s{col.GetParameter(0).IntValue}"
                elif 'ParameterNumber' in par_list: #CPAR
                    comp_name=f"{CPAR_params[col.GetParameter(1).IntValue]}_s{col.GetParameter(0).IntValue}"
                else: #not compensator
                    continue
                #add column for first config only
                if i==1 and not (comp_name==None):
                    df_mc[f"{comp_name}_nominal"] = col.GetParameter(2).DoubleValue #fills same value for all rows, but will be overwritten for rest so doesn't matter
                    df_mc[f"{comp_name}_mean"] = np.nanmean(mcData[:,idx])
                    df_mc[f"{comp_name}_std"] = np.nanstd(mcData[:,idx])
                elif  not (comp_name==None):
                    df_mc.loc[i, f"{comp_name}_nominal"] = col.GetParameter(2).DoubleValue
                    df_mc.loc[i, f"{comp_name}_mean"] = np.nanmean(mcData[:,idx])
                    df_mc.loc[i, f"{comp_name}_std"] = np.nanstd(mcData[:,idx])
    TolReader.Close()
    del zos
    zos=None

its=500
for i in range(1, nconfigs+1):
    run_MC_tol(i, its)

for i in range(1, nconfigs+1):
    get_MC_results(i, df_mc, its)

#make plot of spot radius for each slice with change from nominal
slice_cen*=3600 #convert to arcseconds for plotting
plt.figure('spot size change', figsize=(7,5))
plt.plot(slice_cen, df_mc['RMS_nominal']*1000, 'o', c='tab:cyan', label='Nominal Design')
plt.errorbar(slice_cen, df_mc['RMS_mean']*1000, yerr=df_mc['RMS_std'], fmt='*', c='magenta', label='Toleranced Design')
xlims, ylims = plt.xlim(), plt.ylim() #limits based on results
plt.hlines(y=16.9, xmin=xlims[0], xmax=xlims[1], color='r')#, label='Requirement') 
plt.vlines(x=0, ymin=ylims[0], ymax=ylims[1], lw=1, color='k', alpha=0.2, zorder=0)
plt.xlim(xlims)
#plt.ylim(2,35)
plt.ylim(ylims)
plt.legend(loc='upper center')
plt.xlabel("Slice Position From Center [arcseconds]")
plt.ylabel("RMS Spot Radius [$\\rm \mu m$]")
plt.tight_layout()
#plt.show()
plt.savefig(f"{res_dir}RMS_spot_means_tol_MC_{its}.png", bbox_inches='tight')
plt.close()

def plot_MC_percentile(p, res_dir, iterations):
    plt.figure(f'spot size {p}th', figsize=(7,5))
    plt.plot(slice_cen, df_mc['RMS_nominal']*1000, 'o', c='tab:cyan', label='Nominal Design')
    plt.plot(slice_cen, df_mc[f'RMS_{p}']*1000, '*', c='magenta', label='Toleranced Design')
    xlims, ylims = plt.xlim(), plt.ylim() #limits based on results
    plt.hlines(y=16.9, xmin=xlims[0], xmax=xlims[1], color='r')#, label='Requirement') 
    plt.vlines(x=0, ymin=ylims[0], ymax=ylims[1], lw=1, color='k', alpha=0.2, zorder=0)
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.legend(loc='upper center')
    plt.xlabel("Slice Position From Center [arcseconds]")
    plt.ylabel("RMS Spot Radius [$\\rm \mu m$]")
    plt.tight_layout()
    plt.savefig(f"{res_dir}RMS_spot_p{p}_tol_MC_{iterations}.png", bbox_inches='tight')
    plt.close()

