'''
make_IFU_plots.py
Written by: Maren Cosens
Date: 4/22/26

Description: uses ZOS-API to perform ray tracing on specified slices, wavelengths, and field positions to generate useful plots of image quality across the MIRMOS IFU
-spot diagrams generated for each slice/band with every field position and overlay of three wavelengths spanning band
-RMS spot radius for each band plotted as function of slice position
-may be useful to make cartoon of slice with spot diagrams placed at corresponding field locations
    -colored by wavelength (multiple per band or just center?)
    -perhaps all slices on same figure? may get too busy and want to only include some positions
    (spots fairly consistent across a slice so for now just showing mean per slice - 2D plot could be made following procedure of vignetting plot)

-footprint diagram for all bands with continuous spectra
-plot of spectral and spatial extent on detector as a function of slice position for each band
-plot of wavelength coverage loss as function of field position
-plot of spacing between slices
-plot of vignetting as a function of field position
'''
##import packages
import os, copy
import pandas as pd
import zos_pyclass
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from zemax_functions import get_footprint, get_spot_sizes, make_spot_diagram

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

    
##run functions of interest
if __name__ == '__main__':
    # load local variables
    zos = zos_pyclass.PythonStandaloneApplication()

    ZOSAPI = zos.ZOSAPI
    IFU_System = zos.TheSystem
    #TheApplication = zos.TheApplication
    
    # Setup
    IFU_System.LoadFile(r'C:\Users\mcosens\Documents\Zemax\MIRMOS\IFU\MIRMOS_full_IFS.zos', False)
    #path to save results
    res_dir='C:\\Users\\mcosens\\Documents\\Research_docs\\MIRMOS\\IFU\\'
    
    #set to config1 to start
    IFU_MCE=IFU_System.MCE
    IFU_MCE.SetCurrentConfiguration(1)
    #get system paramaters
    nfields=IFU_System.SystemData.Fields.NumberOfFields
    nwaves=IFU_System.SystemData.Wavelengths.NumberOfWavelengths
    nconfigs=IFU_MCE.NumberOfConfigurations

    ##set shapes and colors for plotting (custom to MIRMOS)
    bands=['i', 'Y', 'J', 'H', 'K']
    band_shapes = ['o', 's', 'p', 'h', 'D']
    band_colors=np.linspace(0.1,1,len(bands))
    #use cycler to color points from sequential colormap
    N = len(bands)*nwaves
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.YlOrRd(np.linspace(0.1,1,N)))

    #get field positions for each slice (should be the same for all configs in a band since only slice position is changing)
    slice_cen = np.full(int(nconfigs/len(bands)), fill_value=np.nan) #center field y position
    field_x = np.full((int(nconfigs/len(bands)), nfields), fill_value=np.nan) #x position of all 7 fields across each slice
    field_y = np.full((int(nconfigs/len(bands)), nfields), fill_value=np.nan) #y position of all 7 fields across each slice
    for c in range(int(nconfigs/len(bands))):
        IFU_MCE.SetCurrentConfiguration(c+1)
        slice_cen[c] = IFU_System.SystemData.Fields.GetField(1).Y
        for f in range(nfields):
            field_x[c,f] = IFU_System.SystemData.Fields.GetField(f+1).X
            field_y[c,f] = IFU_System.SystemData.Fields.GetField(f+1).Y
    slice_cen*=3600 #convert to arcseconds for plotting
    field_x*=3600
    field_y*=3600

    ##spot diagrams (field positions match layout on slice but spacing is not representative)
    '''
    #done already, don't need to repeat
    for c in range(nconfigs):
        band_n=int(c/21)
        if c <21:
            make_spot_diagram(ZOSAPI, IFU_System, c+1, range(1,nfields+1), range(1,nwaves+1), 30, bands[band_n], band_n, res_dir+'\\spot_diagrams\\', norm_waves=False) 
        else:
            make_spot_diagram(ZOSAPI, IFU_System, c+1, range(1,nfields+1), range(1,nwaves+1), 30, bands[band_n], band_n, res_dir+'\\spot_diagrams\\', norm_waves=True)
    
    ##get spot sizes for all fields and wavelengths
    spot_radii = get_spot_sizes(ZOSAPI, IFU_System, range(1, nconfigs+1), range(1,nfields+1), range(1,nwaves+1), 15)[0]
    rms_radii, geo_radii = spot_radii[:,:,:,0], spot_radii[:,:,:,1]

    ##generate plot of spot radii as a function of slice position for each wavelength
    rms_radii[rms_radii==0.0]=np.nan #to filter wavelengths that fall off detector
    
    ##may want to seperate retrieving results and plotting (could save results to fits files in order to close zos application before messing with plots)
    plt.figure(figsize=(7,5))
    for i in range(len(bands)): #loop through each MIRMOS spectral band
        IFU_MCE.SetCurrentConfiguration(1+i*21)
        plt.plot(0,15, 'o', c='none', label=bands[i]+' band') #invisible point to get label in legend
        for j in range(nwaves):
            plt.plot(slice_cen, np.nanmean(rms_radii[i*21:21+i*21,:,j], axis=1), band_shapes[i], label=f"$\\rm{IFU_System.SystemData.Wavelengths.GetWavelength(j+1).Wavelength:.3f} \mu m$") #round wavelength to 3 decimal places
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

    ##detector footprint stats for each band
    #calculate spatial extent on detector for each slice/band
    #use without vignetting to evaluate fraction of wavelength range truncated for each slice
    spec_ycen, spec_xcen = np.full((nconfigs, nfields, nwaves), fill_value=np.nan), np.full((nconfigs, nfields, nwaves), fill_value=np.nan)
    for c in range(nconfigs):
        for f in range(nfields):
            for w in range(nwaves):
                get_footprint(ZOSAPI, IFU_System, config=c+1, field=f+1, wave=w+1, nrays=10, delete_vignetted=False, outpath=res_dir+'footprints\\', outfile=f'footprint_config{c+1}_field{f+1}_wave{w+1}.txt')
    #read in text file and get positions and extent, save to arrays
    for c in range(nconfigs):
        for fd in range(nfields):
            for w in range(nwaves):
                with open(f"{res_dir}\\footprints\\footprint_config{c+1}_field{fd+1}_wave{w+1}.txt", 'r', encoding='utf-16') as f:
                    lines = f.readlines()
                    for line in lines:
                        #get lines with '=' and split on that
                        if '=' in(line):
                            param, value = line.split('=')
                            if param.strip()=='Ray X Center':
                                spec_xcen[c,fd,w] = float(value.strip())
                            elif param.strip()=='Ray Y Center':
                                spec_ycen[c,fd,w] = float(value.strip())
    #use positions to plot footprint diagrams for each band, color by slice (fill is approximate)
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20b(np.concatenate([np.linspace(0.5,1,11), np.linspace(0.5,0,11)])))
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,10), layout='tight')
    col_num=np.linspace(0.1,1,21)
    for i in range(len(bands)-1):
        row = i%2
        col = int(i/2)
        ax[row, col].set_title(f"{bands[i+1]} band")
        for k in range((i+1)*21,(i+1)*21+21):
            ax[row, col].fill_between(x=[np.nanmin(spec_xcen[k,:,:]), np.nanmax(spec_xcen[k,:,:])], y1=np.nanmin(spec_ycen[k,:,:]), y2=np.nanmax(spec_ycen[k,:,:]))
        ax[row, col].set_xlim(-18.4, 18.4)
        ax[row, col].set_ylim(-18.4, 18.4)
        ax[row, col].set_xlabel("X Position on Detector [mm]")
        ax[row, col].set_ylabel("Y Position on Detector [mm]")
    plt.tight_layout()
    plt.savefig(f"{res_dir}footprint_diagrams.png", bbox_inches='tight')
    plt.close()

    #get spatial and spectral extent from arrays and make plot of values by band (later can make 2D plot of spatial vs spectral extent for each slice)
    spec_extent = np.sqrt((spec_ycen[:,1,2]-spec_ycen[:,1,0])**2 + (spec_xcen[:,1,2]-spec_xcen[:,1,0])**2) #use central field
    spatial_extent = np.sqrt((spec_ycen[:,-1,2]-spec_ycen[:,-2,2])**2 + (spec_xcen[:,-1,2]-spec_xcen[:,-2,2])**2) #use central wavelength
    #later add plot of space between slices by band
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.YlOrRd(np.linspace(0.1,1,len(bands))))
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7,10), layout='tight', sharex=True)
    for i in range(len(bands)): #loop through each MIRMOS spectral band
        ax[0].plot(slice_cen, spatial_extent[i*21:21+i*21], band_shapes[i], label=bands[i]+' band')
        ax[1].plot(slice_cen, spec_extent[i*21:21+i*21], band_shapes[i], label=bands[i]+' band')
    ax[1].set_ylim(26,30) #leave off i-band since spectral extent is not meaningful
    ax[0].set_ylabel("Spatial Extent on Detector [mm]")
    secax = ax[0].secondary_yaxis('right', functions=(lambda x: (x - 1.222) / 1.222 * 100,
                                                      lambda y: y / 100 * 1.222 + 1.222))
    xlims = ax[0].get_xlim() #limits based on results
    ax[0].hlines(y=1.222, xmin=xlims[0], xmax=xlims[1], color='gray', alpha=0.7, ls='--', zorder=0)
    ax[0].set_xlim(xlims)
    secax.set_ylabel("Spatial Extent Deviation from Nominal [%]")
    ax[1].set_ylabel("Spectral Extent on Detector [mm]")
    ax[1].set_xlabel("Slice Position From Center [arcseconds]")
    ax[1].set_ylabel("Spectral Extent on Detector [mm]")
    ax[0].legend(ncols=5, loc='upper center')
    plt.tight_layout()
    plt.savefig(f"{res_dir}detector_extent.png", bbox_inches='tight')
    plt.close()

    #space between slices for each band
    #start with plotting central wavelength, then check it's the same at edge wavelengths
    spec_gap = np.full((nconfigs, nwaves), fill_value=np.nan) #need one less config since this is difference between adjacent slices
    for b in range(len(bands)):
        for i in range(1,len(slice_cen)):
            if slice_cen[i]>0: 
                #previous slice field6[5] - current slice field7[6]
                spec_gap[(i-1)+b*21,:] = spec_xcen[i-1+b*21,5,:]-spec_xcen[i+b*21,6,:] - (rms_radii[i-1+b*21,5,:]+rms_radii[i+b*21,6,:])/1000 #subtract spot radius to get gap between slices rather than centroids
            elif i==11:
                spec_gap [(i-1)+b*21,:] = spec_xcen[i+b*21,5,:]-spec_xcen[b*21,6,:] - (rms_radii[i+b*21,5,:]+rms_radii[b*21,6,:])/1000#first negative slice needs to be compared to center slice
            else:
                #current slice field6[5] - previous slice field7[6]
                spec_gap[(i-1)+b*21,:] = spec_xcen[i+b*21,5,:]-spec_xcen[i-1+b*21,6,:] - (rms_radii[i+b*21,5,:]+rms_radii[i-1+b*21,6,:])/1000 #subtract spot radius to get gap between slices rather than centroids
    spec_gap = spec_gap[~np.isnan(spec_gap).any(axis=1)] #drop nan values from spec_gap (no gap for center slice)

    #plot results
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.YlOrRd(np.linspace(0.1,1,N)))
    avg_cen = np.concatenate([np.linspace(0.5*.84,9.5*.84, int((len(slice_cen)-1)/2)), np.linspace(-0.5*.84,-9.5*.84, int((len(slice_cen)-1)/2))])
    plt.figure(figsize=(7,5))
    for i in range(len(bands)): #loop through each MIRMOS spectral band
        IFU_MCE.SetCurrentConfiguration(1+i*21)
        plt.plot(0,15, 'o', c='none', label=bands[i]+' band') #invisible point to get label in legend
        for j in range(nwaves):
            if i==0: #imaging chanel has different pixle scale
                plt.plot(avg_cen, spec_gap[i*20:20+i*20,j]/0.015, band_shapes[i], label=f"$\\rm{IFU_System.SystemData.Wavelengths.GetWavelength(j+1).Wavelength:.3f} \mu m$") #round wavelength to 3 decimal places
            else:
                plt.plot(avg_cen, spec_gap[i*20:20+i*20,j]/0.018, band_shapes[i], label=f"$\\rm{IFU_System.SystemData.Wavelengths.GetWavelength(j+1).Wavelength:.3f} \mu m$") #round wavelength to 3 decimal places
    plt.xlabel("Slice Position From Center [arcseconds]")
    plt.ylabel("Gap Between Slices [pixels]")
    plt.legend(ncols=5, loc='upper center')
    plt.tight_layout()
    plt.savefig(f"{res_dir}detector_spacing.png", bbox_inches='tight')
    plt.close()

    ##wavelength truncation by band/slice
    #get distance above or below edge of detector and compare to spectral extent of that slice to get fraction truncated and what wavelengths are covered
    spectral_coverage = np.zeros(nconfigs) #negative=red end truncated, positive=blue end truncated, 0=no truncation, abs value is fraction of wavelengths truncated based on distance from edge of detector and spectral extent of slice
    for i in range(nconfigs):
        if spec_ycen[i,1,0]>18.4: #blue end truncated
            spectral_coverage[i] = (spec_ycen[i,1,0]-18.4)/spec_extent[i] #fraction of wavelengths truncated based on distance from edge of detector and spectral extent of slice
        elif spec_ycen[i,1,2]<-18.4: #red end truncated
            spectral_coverage[i] = (18.4+spec_ycen[i,1,2])/spec_extent[i] #fraction of wavelengths truncated based on distance from edge of detector and spectral extent of slice
    #plot results
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.YlOrRd(np.linspace(0.325,1,len(bands)-1)))
    plt.figure(figsize=(7,5)) 
    for i in range(1,len(bands)): #loop through each MIRMOS spectral bands
        plt.plot(slice_cen, spectral_coverage[i*21:21+i*21]*100, band_shapes[i], label=bands[i]+' band')
    plt.xlabel("Slice Position From Center [arcseconds]")
    plt.ylabel("% of Wavelength Range Truncated")
    plt.legend(ncols=5, loc='upper center')
    plt.tight_layout()
    plt.savefig(f"{res_dir}spectral_truncation.png", bbox_inches='tight')
    plt.close()

    #generate table of spectral coverage by slice and band from spectral_coverge array
    spectral_table = pd.DataFrame(spectral_coverage.reshape(21, 5, order='F'), columns=[f"{bands[i]}_trunc_frac" for i in range(len(bands))])
    for i in range(1,len(bands)):
        IFU_MCE.SetCurrentConfiguration(1+i*21)
        band_start = IFU_System.SystemData.Wavelengths.GetWavelength(1).Wavelength
        band_end = IFU_System.SystemData.Wavelengths.GetWavelength(nwaves).Wavelength
        spectral_table[f"{bands[i]}_start"] = band_start
        spectral_table[f"{bands[i]}_end"] = band_end
        for j in range(len(spectral_table)):
            if spectral_table[f"{bands[i]}_trunc_frac"][j]<0: #red end truncated
                spectral_table.loc[j, f"{bands[i]}_end"] = band_end + spectral_table[f"{bands[i]}_trunc_frac"][j]*(band_end - band_start)
            elif spectral_table[f"{bands[i]}_trunc_frac"][j]>0: #blue end truncated
                spectral_table.loc[j, f"{bands[i]}_start"] = band_start + spectral_table[f"{bands[i]}_trunc_frac"][j]*(band_end - band_start)
    spectral_table.to_csv(f"{res_dir}spectral_coverage_table.csv")
    
    ##Evaluate footprint diagrams at pupil for each band to get centration and radius
    for c in range(nconfigs):
        get_footprint(ZOSAPI, IFU_System, config=c+1, field=0, wave=0, nrays=10, surface=66, delete_vignetted=False, outpath=res_dir+'footprints\\', outfile=f'footprint_config{c+1}_allfields_allwaves_pupil.txt')
    #read in text file and get values for center and radii
    pupil_xcen = np.full(nconfigs, fill_value=np.nan)
    pupil_ycen = np.full(nconfigs, fill_value=np.nan)
    pupil_xrad = np.full(nconfigs, fill_value=np.nan)
    pupil_yrad = np.full(nconfigs, fill_value=np.nan)
    for c in range(nconfigs):
        with open(f"{res_dir}\\footprints\\footprint_config{c+1}_allfields_allwaves_pupil.txt", 'r', encoding='utf-16') as f:
            lines = f.readlines()
            for line in lines:
                if '=' in(line):
                    param, value = line.split('=')
                    if param.strip()=='Ray X Center':
                        pupil_xcen[c] = float(value.strip())
                    elif param.strip()=='Ray Y Center':
                        pupil_ycen[c] = float(value.strip())
                    elif param.strip()=='Ray X Half Width':
                        pupil_xrad[c] = float(value.strip())
                    elif param.strip()=='Ray Y Half Width':
                        pupil_yrad[c] = float(value.strip())
    #generate plot of pupil centration and radius by slice position and band
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(7,10), layout='tight', sharex=True)
    for i in range(len(bands)): #loop through each MIRMOS spectral band
        ax[0].plot(slice_cen, np.abs(pupil_xcen[i*21:21+i*21]), band_shapes[i], mfc='none', mec=plt.cm.YlOrRd(band_colors[i]))#, label=bands[i]+' band') 
        ax[0].plot(slice_cen, np.abs(pupil_ycen[i*21:21+i*21]), band_shapes[i], c=plt.cm.YlOrRd(band_colors[i]), label=bands[i]+' band')
        ax[1].plot(slice_cen, pupil_xrad[i*21:21+i*21], band_shapes[i], mfc='none', mec=plt.cm.YlOrRd(band_colors[i]))#label=bands[i]+' band')
        ax[1].plot(slice_cen, pupil_yrad[i*21:21+i*21], band_shapes[i], c= plt.cm.YlOrRd(band_colors[i]), label=bands[i]+' band') 
    ax[0].set_ylabel("Pupil Decenter [mm]")
    ax[1].set_ylabel("Pupil Radius [mm]")
    ax[1].set_xlabel("Slice Position From Center [arcseconds]")
    xlims, ylims = ax[1].get_xlim(), ax[1].get_ylim() #limits based on results
    ax[1].hlines(y=58.5, xmin=xlims[0], xmax=xlims[1], color='r', alpha=0.7, ls='solid', zorder=0, label='Lyot stop') #nominal Lyot stop radius in K radius
    ax[1].plot(xlims[0]-10,ylims[0]-10, 'o',mfc='none', mec='k', label='x direction')
    ax[1].plot(xlims[0]-10,ylims[0]-10, 'o', c='k', label='y direction')
    ax[1].set_xlim(xlims)
    ax[1].set_ylim(ylims)
    #customize order of labels in legend
    handles, labels = ax[1].get_legend_handles_labels()
    order = [0, 6, 1, 7, 2, 5, 3, 4]
    ax[0].legend([handles[i] for i in order], [labels[i] for i in order], ncols=5, loc='upper center')
    plt.tight_layout()
    plt.savefig(f"{res_dir}pupil_metrics.png", bbox_inches='tight')
    plt.close()
    '''

    ##plot of vignetting as function of field position
    #need to add Lyot stop aperture for accurate vignetting
    #will need to repeat footprint analysis with delete_vignetted=True (use only central wavelength to remove impact of spectra falling off detector)
    vignet_frac = np.full((nconfigs, nfields), fill_value=np.nan)
    for c in range(nconfigs):
        for f in range(nfields):
            get_footprint(ZOSAPI, IFU_System, config=c+1, field=f+1, wave=2, nrays=10, delete_vignetted=True, outpath=res_dir+'footprints\\', outfile=f'footprint_config{c+1}_field{f+1}_wave2_vignetting.txt')
    #read in text file and get percent of rays through
    for c in range(nconfigs):
        for fd in range(nfields):
                with open(f"{res_dir}\\footprints\\footprint_config{c+1}_field{fd+1}_wave2_vignetting.txt", 'r', encoding='utf-16') as f:
                    lines = f.readlines()
                    line = lines[-1]
                    param, val = line.split('=')
                    val_num, val_unit = val.split('%')
                    vignet_frac[c,fd] = 1- float(val_num.strip())/100 #convert to fraction of rays vignetted rather than through
    #plot as function of x,y field position
    vign_nans=copy.copy(vignet_frac)
    vign_nans[vign_nans==0]=np.nan #to filter out points with no vignetting for plotting
    max_vign=np.nanmax(vignet_frac)
    #interpolate between field points
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(10,10))
    for i in range(len(bands)):
        if i==0:
            row, col = 2,0 #put at the bottom since this is the least relevant for vignetting
        else:
            j=i-1
            col = j%2
            row = int(j/2)
        ax[row, col].set_title(f"{bands[i]} band")
        ax[row,col].pcolormesh(field_x, field_y, vign_nans[i*21:21+i*21,:], cmap='bone_r', vmin=0, vmax=max_vign, shading = 'gouraud') #not the smoothest but seems to be best option without dealing with field positions not being sequential or regular grid
        ax[row,col].set_xlabel("Field X [arcseconds]")
        ax[row,col].set_ylabel("Field Y [arcseconds]")
    ax[2, 1].axis('off')
    sm = plt.cm.ScalarMappable(cmap='bone_r', norm=plt.Normalize(vmin=0.01, vmax=max_vign))
    sm.set_array([])
    fig.colorbar(sm, ax=ax.ravel().tolist(), label='Fraction of Rays Vignetted', pad=0.05, location='right')
    fig.subplots_adjust(wspace=0.25, hspace=0.3, top=0.95, bottom=0.05, left=0.1, right=0.75)
    plt.savefig(f"{res_dir}vignetting_smooth.png", bbox_inches='tight')
    plt.close()
    #discrete field points
    fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(10,10))
    for i in range(len(bands)):
        if i==0:
            row, col = 2,0 #put at the bottom since this is the least relevant for vignetting
        else:
            j=i-1
            col = j%2
            row = int(j/2)
        ax[row, col].set_title(f"{bands[i]} band")
        ax[row, col].scatter(field_x.flatten(), field_y.flatten(), c=vignet_frac[i*21:21+i*21,:].flatten(), marker='s', cmap='bone_r', vmin=0.01, vmax=max_vign, alpha=0.8)
        ax[row,col].set_xlabel("Field X Position [arcseconds]")
        ax[row,col].set_ylabel("Field Y Position [arcseconds]")  
    ax[2, 1].axis('off')
    fig.colorbar(sm, ax=ax.ravel().tolist(), label='Fraction of Rays Vignetted', pad=0.05, location='right')
    fig.subplots_adjust(wspace=0.25, hspace=0.3, top=0.95, bottom=0.05, left=0.1, right=0.75)
    plt.savefig(f"{res_dir}vignetting_discrete.png", bbox_inches='tight')
    plt.close()
    ##single band vignetting diagram
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6.5,5))
    i=4
    ax.set_title(f"{bands[i]}-band Vignetting")
    ax.scatter(field_x.flatten(), field_y.flatten(), c=vignet_frac[i*21:21+i*21,:].flatten()*100, marker='s', cmap='bone_r', vmin=0.01*100, vmax=max_vign*100, alpha=1)
    ax.set_xlabel("Field X Position [arcseconds]")
    ax.set_ylabel("Field Y Position [arcseconds]")
    ax.set_xlim(np.nanmin(field_x)-.5, np.nanmax(field_x)+.5)
    ax.set_ylim(np.nanmin(field_y)-.5, np.nanmax(field_y)+.5)
    sm = plt.cm.ScalarMappable(cmap='bone_r', norm=plt.Normalize(vmin=0.01*100, vmax=max_vign*100))
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='% of Rays Vignetted', pad=0.05, location='right') 
    plt.tight_layout()
    plt.savefig(f"{res_dir}vignetting_{bands[i]}band.png", bbox_inches='tight')
    plt.close()
    #make 1-D cuts?

    # close server instance of OpticStudio
    del zos
    zos = None


    