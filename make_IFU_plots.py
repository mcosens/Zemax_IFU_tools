'''
make_IFU_plots.py
Written by: Maren Cosens
Date: 4/22/26

Description: uses ZOS-API to perform ray tracing on specified slices, wavelengths, and field positions to generate useful plots of image quality across the MIRMOS IFU
-may be useful to make cartoon of slice with spot diagrams placed at corresponding field locations
    -colored by wavelength (multiple per band or just center?)
    -perhaps all slices on same figure? may get too busy and want to only include some positions
-plot of vignetting as a function of field position
-plot of wavelength coverage loss as function of field position
-footprint plots for all bands (very slow to make in zemax with this many configurations)
'''
##import packages
import zos_pyclass
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from System import Enum, Int32, Double

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

    
##functions for various plots
def make_spot_diagram(system, config, fields, waves, nrays, band, band_n, outpath=".\\", norm_waves=True):
    """
    Perform a ray trace using the ZOS-API and generate spot diagrams for a given configurations, fields, and wavelengths.

    Parameters:
    - system: The ZOS system object
    - config: Configuration number
    - fields: List of field numbers
    - waves: List of wavelength numbers
    - nrays: Number of rays for the analysis
    - band: String name of band
    - band_n: band number (to derive slice number)
    - outpath: filepatch to save plot (default=".\\")
    - norm_waves: bool; default = True; indicates whether to shift spot diagram centroids to 0,0

    Returns:
    - spot diagrams saved to outpath+'spot_diagram_slice{slice}{band}.png'
    """
    system.MCE.SetCurrentConfiguration(config) #switch to desired configuration
    slice=config-21*band_n
    #ray tracing adapted from 'PythonStandalone_22_seq_spot_diagram.py'
    # Set up Batch Ray Trace
    raytrace = system.Tools.OpenBatchRayTrace()
    nsur = system.LDE.NumberOfSurfaces
    max_rays = nrays
    normUnPolData = raytrace.CreateNormUnpol((max_rays + 1) * (max_rays + 1), ZOSAPI.Tools.RayTrace.RaysType.Real, nsur)

    # Define batch ray trace constants
    max_wave = len(waves)
    num_fields=len(fields)
    field_x_ar = [system.SystemData.Fields.GetField(int(f)).X for f in fields]
    field_y_ar = [system.SystemData.Fields.GetField(int(f)).Y for f in fields]
    hy_ar = field_y_ar/(np.abs(np.max(field_y_ar))) #get normalization of field points
    hx_ar = field_x_ar/(np.abs(np.max(field_x_ar)))
    
    # Initialize x/y image plane arrays
    x_ar = np.empty((num_fields, max_wave, ((max_rays + 1) * (max_rays + 1))))
    y_ar = np.empty((num_fields, max_wave, ((max_rays + 1) * (max_rays + 1))))

    # Determine maximum field in X,Y
    max_field_y, max_field_x =np.max(field_y_ar), np.max(field_x_ar)

    if system.SystemData.Fields.GetFieldType() == ZOSAPI.SystemData.FieldType.Angle:
        field_type = 'Angle'
    elif system.SystemData.Fields.GetFieldType() == ZOSAPI.SystemData.FieldType.ObjectHeight:
        field_type = 'Height'
    elif system.SystemData.Fields.GetFieldType() == ZOSAPI.SystemData.FieldType.ParaxialImageHeight:
        field_type = 'Height'
    elif system.SystemData.Fields.GetFieldType() == ZOSAPI.SystemData.FieldType.RealImageHeight:
        field_type = 'Height'
    
    #set-up figure prior to ray trace for each field
    rows, cols = len(set(hy_ar)), len(set(hx_ar))
    colors = ('y', 'm', 'c')
    fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=(cols*2,rows*2), layout='tight') #set up subplots so each field is square, but doesn't reproduce aspect ratio of slice...(work on this later)
    drop_ax=np.full((rows, cols), fill_value=True) #where this remains True after loop, drop these axes
    for field in fields:
        #get row and column based on field position relative to others
        hy, hx = hy_ar[field-1], hx_ar[field-1]
        row = np.argwhere(sorted(set(hy_ar))==hy)[0][0]
        col = np.argwhere(sorted(set(hx_ar))==hx)[0][0]
        drop_ax[row, col] = False
        ax[row, col].set_title(f'{system.SystemData.Fields.GetField(field).Comment}')
        for wave in waves:
            # Adding Rays to Batch, varying normalised object height hy
            normUnPolData.ClearData()
            #for i = 1:((max_rays + 1) * (max_rays + 1))
            for i in range(1, (max_rays + 1) * (max_rays + 1) + 1):
                px = np.random.random() * 2 - 1
                py = np.random.random() * 2 - 1
                while (px*px + py*py > 1):
                    py = np.random.random() * 2 - 1
                normUnPolData.AddRay(wave, hx_ar[field -1], hy_ar[field -1], px, py, Enum.Parse(ZOSAPI.Tools.RayTrace.OPDMode, "None"))
            
            raytrace.RunAndWaitForCompletion()
            # Read batch raytrace and display results
            normUnPolData.StartReadingResults()
            
            # Python NET requires all arguments to be passed in as reference, so need to have placeholders
            sysInt = Int32(1)
            sysDbl = Double(1.0)
            
            output = normUnPolData.ReadNextResult(sysInt, sysInt, sysInt,
                            sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl);

            while output[0]:                                                    # success
                if ((output[2] == 0) and (output[3] == 0)):                     # ErrorCode & vignetteCode
                    x_ar[field - 1, wave - 1, output[1] - 1] = output[4]   # X
                    y_ar[field - 1, wave - 1, output[1] - 1] = output[5]   # Y
                output = normUnPolData.ReadNextResult(sysInt, sysInt, sysInt,
                            sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl, sysDbl);
            
            ##filter results for plots
            #remove zeros in results (throws off plot scale and outlier filtering)
            x_ar[field-1, wave-1, :][x_ar[field-1, wave-1, :]==0], y_ar[field-1, wave-1, :][y_ar[field-1, wave-1, :]==0] = np.nan, np.nan
            x_med = np.nanmedian(x_ar[field-1, wave-1, :])
            y_med = np.nanmedian(y_ar[field-1, wave-1, :])
            x_std = np.nanstd(x_ar[field-1, wave-1, :])
            y_std = np.nanstd(y_ar[field-1, wave-1, :])
            #remove significant outliers in results (throws off plot scale)
            y_ar[field-1, wave-1, :][y_ar[field-1, wave-1, :] > y_med+5*y_std] = np.nan
            y_ar[field-1, wave-1, :][y_ar[field-1, wave-1, :] < y_med-5*y_std] = np.nan
            x_ar[field-1, wave-1, :][x_ar[field-1, wave-1, :] > x_med+5*x_std] = np.nan
            x_ar[field-1, wave-1, :][x_ar[field-1, wave-1, :] < x_med-5*x_std] = np.nan
            if norm_waves:
                #apply shift
                x_ar[field-1, wave-1, :] -= x_med
                y_ar[field-1, wave-1, :] -= y_med
            #create plot item
            temp = ax[row,col].plot(np.squeeze(x_ar[field - 1, wave - 1, :]), np.squeeze(y_ar[field - 1, wave - 1, :]), '.', ms = 1, color = colors[wave - 1], label=f"$\\rm{system.SystemData.Wavelengths.GetWavelength(wave).Wavelength:.3f} \mu m$") 
        #add circle for requirement (at centroid)
        cx, cy = np.nanmedian(x_ar[field - 1, :, :]), np.nanmedian(y_ar[field - 1, :, :]) 
        circle = plt.Circle((cx,cy), (2.355/2)*0.0169, fill=False, edgecolor='red', linewidth=1, linestyle='--')
        ax[row, col].add_patch(circle)
    #drop axes not used
    bad_ax=np.argwhere(drop_ax)
    for i in range(bad_ax.shape[0]):
        ax[bad_ax[i][0], bad_ax[i][1]].remove()
    plt.suptitle('Spot Diagram, Slice: %s' % (slice))
    #make legend with only unique labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), markerscale=5, borderaxespad=5)
    #plt.subplots_adjust(wspace=0.8)
    plt.draw()
    plt.savefig(f"{outpath}spot_diagram_slice{slice}{band}.png", bbox_inches='tight')
    plt.close()
    system.Tools.CurrentTool.Close() #need to add this in order to do a subsequent ray trace

def get_spot_sizes(system, configs, fields, waves, nrays):
    """
    Computes spot sizes (RMS and geometric) for given configurations, fields, and wavelengths using ZOS-API spot analysis.

    Parameters:
    - system: The ZOS system object
    - configs: List of configuration numbers
    - fields: List of field numbers
    - waves: List of wavelength numbers
    - nrays: Number of rays for the analysis

    Returns:
    - results: numpy array of shape (len(configs), len(fields), len(waves), 2) containing RMS and Geo spot sizes
    """
    #IFU_MCE=IFU_System.MCE
    spot = system.Analyses.New_Analysis_SettingsFirst(ZOSAPI.Analysis.AnalysisIDM.StandardSpot)
    spot_setting = spot.GetSettings().__implementation__ #.__implementation__ added to deal with errors from using PythonNet3
    spot_setting.Surface.UseImageSurface()
    spot_setting.ReferTo = ZOSAPI.Analysis.Settings.Spot.Reference.Centroid
    spot_setting.Patterns = ZOSAPI.Analysis.Settings.Spot.Patterns.Dithered 
    spot_setting.RayDensity = nrays
    results = np.full((len(configs), len(fields), len(waves), 2), fill_value=np.nan) #empty array for RMS and Geo radius for each field/wavelength
    slice_center = np.full(len(configs), fill_value=np.nan)
    for k in range(len(configs)):
        c=configs[k]
        system.MCE.SetCurrentConfiguration(c)
        slice_center[k] = system.SystemData.Fields.GetField(1).Y
        for i in range(len(fields)):
            for j in range(len(waves)):
                f, w = fields[i], waves[j]
                spot_setting.Field.SetFieldNumber(f)
                spot_setting.Wavelength.SetWavelengthNumber(w)
                spot.ApplyAndWaitForCompletion()
                spot_results = spot.GetResults()
                results[k, i, j, :] = spot_results.SpotData.GetRMSSpotSizeFor(f, w), spot_results.SpotData.GetGeoSpotSizeFor(f, w)
    return results, slice_center

##run functions of interest
if __name__ == '__main__':
    # load local variables
    zos = zos_pyclass.PythonStandaloneApplication()

    ZOSAPI = zos.ZOSAPI
    IFU_System = zos.TheSystem
    TheApplication = zos.TheApplication
    
    # Setup
    IFU_System.LoadFile(r'C:\Users\mcosens\Documents\Zemax\MIRMOS\IFU\MIRMOS_full_IFS.zos', False)
    #path to save results
    res_dir='C:\\Users\\mcosens\\Documents\\Research_docs\\MIRMOS\\IFU\\'
    
    #set to config1 for test
    IFU_MCE=IFU_System.MCE
    IFU_MCE.SetCurrentConfiguration(1)
    #get system paramaters
    nfields=IFU_System.SystemData.Fields.NumberOfFields
    nwaves=IFU_System.SystemData.Wavelengths.NumberOfWavelengths
    nconfigs=IFU_MCE.NumberOfConfigurations

    ##set shapes and colors for plotting (custom to MIRMOS)
    band_shapes = ['o', 's', 'p', 'h', 'D']
    bands=['i', 'Y', 'J', 'H', 'K']
    #use cycler to color points from sequential colormap
    N = len(bands)*nwaves
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.YlOrRd(np.linspace(0.1,1,N)))

    ##spot diagrams (field positions match layout on slice but spacing is not representative)
    '''
    #done already, don't need to repeat
    for c in range(nconfigs):
        band_n=int(c/21)
        if c <21:
            make_spot_diagram(IFU_System, c+1, range(1,nfields+1), range(1,nwaves+1), 30, bands[band_n], band_n, res_dir+'\\spot_diagrams\\', norm_waves=False) 
        else:
            make_spot_diagram(IFU_System, c+1, range(1,nfields+1), range(1,nwaves+1), 30, bands[band_n], band_n, res_dir+'\\spot_diagrams\\', norm_waves=True)
    
    ##get spot sizes for all fields and wavelengths
    spot_radii, slice_cen = get_spot_sizes(IFU_System, range(1, nconfigs+1), range(1,nfields+1), range(1,nwaves+1), 15)
    rms_radii, geo_radii = spot_radii[:,:,:,0], spot_radii[:,:,:,1]

    ##generate plot of spot radii as a function of slice position for each wavelength
    rms_radii[rms_radii==0.0]=np.nan #to filter wavelengths that fall off detector
    
    ##may want to seperate retrieving results and plotting (could save results to fits files in order to close zos application before messing with plots)
    plt.figure(figsize=(7,5))
    for i in range(5): #loop through each MIRMOS spectral band
        IFU_MCE.SetCurrentConfiguration(1+i*21)
        plt.plot(0,15, 'o', c='none', label=bands[i]+' band') #invisible point to get label in legend
        for j in range(nwaves):
            plt.plot(slice_cen[:21]*3600, np.nanmean(rms_radii[i*21:21+i*21,:,j], axis=1), band_shapes[i], label=f"$\\rm{IFU_System.SystemData.Wavelengths.GetWavelength(j+1).Wavelength:.3f} \mu m$") #round wavelength to 3 decimal places
    xlims, ylims = plt.xlim(), plt.ylim() #limits based on results
    plt.hlines(y=16.9, xmin=xlims[0], xmax=xlims[1], color='r')#, label='Requirement') 
    plt.vlines(x=0, ymin=ylims[0], ymax=ylims[1], lw=1, color='k', alpha=0.2, zorder=0)
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.legend(ncols=5, loc='upper center')
    plt.xlabel("Slice Position From Center [arcseconds]")
    plt.ylabel("RMS Spot Radius [$\\rm \mu m$]")
    plt.tight_layout()
    plt.savefig(f"{res_dir}RMS_spot_means.png", bbox_inches='tight')
    plt.close()
    '''
    ##detector footprint diagrams for each band
    #start with matching Zemax to verify, then make spectra continuous in wavelength
    #use to evaluate fraction of wavelength range truncated for each slice
    IFU_System.Analyses.New_Analysis(ZOSAPI.Analysis.AnalysisIDM.FootprintSettings) #Footprint Diagram not yet built in
    #calculate spatial extent on detector for each slice/band
    #space between slices for each band


    ##footprint diagrams at pupil for each band
    #use to calculate vignetting as function of field position

    ##footprints at pupil mirrors (use to set aperture for each by modifying MCE)

    # close server instance of OpticStudio
    del zos
    zos = None