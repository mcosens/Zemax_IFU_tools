# Zemax_IFU_tools
scripts to interface with the ZOS-API and automate aspects of the design and evaluation of the MIRMOS IFU

**WARNING:** These scripts are currently a work in-progress and currently are quite specific to the design of the MIRMOs IFU. Generalizing and expanding their functionality is TBD

## Description
**Main Scripts**
1. `make_full_model.py`: uses full spectrograph zos file and single band IFU model to create new combined model with all spatial slices and spectral bands
2. `make_IFU_plots.py`: perform ray tracing and generate spot / footprint diagrams using full system model to evalueate performance of the IFU at all bands / slices
3. `set_pupil_mirror_aperture.py`: perform footprint analysis for each pupil mirror across all bands to determine required aperture to provide 5% CA
4. `tolerance_IFU.py`: perform sensitivity analysis on each slice to set tolerance parameters prior to running full Monte Carlo analysis to evaluate expected performance of full IFU

**Other Functions**

1. `zos_pyclass.py`: class definition allowing interface with the ZOS-API in "standalone mode" (adapted ZOS-API example scripts)
2. `zemax_functions.py`: container for various functions used elsewhere (e.g., to make footprint diagrams, convert from local to global coordinates, etc.)

## Citing
If you find these tools useful and make use of them in your work, please cite the SPIE proceedings Cosens et al. (2026)
*(link to be added when published to the SPIE Digital Library)*