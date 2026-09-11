# SEPIA (SuscEptibility mapping PIpeline tool for phAse images)

![sepia logo](https://sepia-documentation.readthedocs.io/en/latest/_static/logo.png)

## Introduction  

**SEPIA** is a tool providing a graphical user interface to build data processing pipeline of quantitative susceptibility mapping (QSM) in Matlab.

The current GUI version is built to access the following toolboxes:
- [MEDI (updated Jan 15, 2020)](http://weill.cornell.edu/mri/pages/qsm.html), 
- [STI Suite (v3.0)](https://chunleiliulab.github.io/software.html),
- [FANSI (v3.0, released on 2021.10.15, i.e., commit b6ac1c9e)](https://gitlab.com/cmilovic/FANSI-toolbox/-/tree/b6ac1c9ea03380722ebe25a6dbef33fff4ea3700),  
- [SEGUE](https://xip.uclb.com/i/software/SEGUE.html), and 
- [nonlinear dipole inversion (NDI)](https://github.com/polakd/NDI_Toolbox),
- [mritools (ROMEO/CLEARSWI) (v4.6.1)](https://github.com/korbinian90/CompileMRI.jl/releases/tag/v4.6.1),
- [MRI Susceptibility Calculation Methods, accessed 12 September 2019](https://xip.uclb.com/product/mri_qsm_tkd),
- HEIDI (auto-download script available - see the documentation), and
- [Chi-separation toolbox](https://github.com/SNU-LIST/chi-separation).

SEPIA provides two key features for QSM processing:  
1. mix-and-match methods from different toolboxes to build your own QSM processing pipeline,
2. graphical user interface to easily adjust parameters of different algorithms.

SEPIA is designed to provide a platform for easy access to different QSM processing methods in the field. To achieve this, most of the codes were written for data flow and algorithm parameter control. Through SEPIA, we hope researchers who are not expert in QSM will also be able to use QSM for their research.

**For better readability, the documentation of SEPIA has moved to https://sepia-documentation.readthedocs.io/.**  

## Terms of use
All the original codes and methods developed for **SEPIA** are under MIT license. You can check [the license file](https://github.com/kschan0214/Sepia/blob/master/LICENSE) for more information. For the terms of use of the toolboxes related to this work, their own license applied and please check the corresponding license file(s) in each toolbox for more information. 

If you use SEPIA in your research, please cite the following article:

[Chan, K.-S., Marques, J.P., 2021. SEPIA—Susceptibility mapping pipeline tool for phase images. Neuroimage 227, 117611.](https://doi.org/10.1016/j.neuroimage.2020.117611)  

As well as any related papers in your processing pipeline. 

If you encounter a bug in SEPIA, please report to [github page](https://github.com/kschan0214/sepia/issues). 

If you have a more general question regarding the usage of SEPIA and/or other QSM questions, please make use of [github page](https://github.com/kschan0214/sepia/discussions).


## Update notes  

For full update log, please visit https://sepia-documentation.readthedocs.io/en/latest/getting_started/Release-note.html.

### 1.3.0 (current, commit 3ad47fd)

Thank you to everyone who contributed to this major release! Special thanks to Patrick, whose engagement and input made many of this update's features possible.

> **Upgrade notes**
> `SpecifyToolboxesDirectory.m` and `SpecifyAtlasDirectory.m` are no longer tracked in git (see "Housekeeping" below). No action needed - your existing local copies are untouched - but `git status` will now show them as untracked instead of unmodified.

**New QSM methods & toolboxes**
* Added support for the χ-separation (Chi-separation) toolbox as a new QSM add-on (paramagnetic/diamagnetic susceptibility separation via Chi-sepnet, chi_sep_MEDI and chi_sep_iLSQR; requires ONNX checkpoint files and the Deep Learning Toolbox Converter for ONNX Model Format support package) Thanks to Taechang and the SNU team for making this possible!
* Added HEIDI as a dipole inversion method, available both as the two-stage "LSQR+HEIDI" pipeline and as a "Streaking reduction by HEIDI" post-processing option that can be applied on top of any other QSM dipole-inversion method's output. Thanks to Fahad and Ferdinand for sharing the code and libraries that made this possible!
* Updated the `mu2` parameter handling for FANSI. Thanks to Sebastian for the bug fix.
* New `download_FANSI_toolbox.m` script to automatically download a pinned FANSI-toolbox commit and register it in `SpecifyToolboxesDirectory.m`
* `HEIDI_HOME` and `ChiSepNet_HOME` are now configured centrally in `SpecifyToolboxesDirectory.m` (editable via the Utility tab's Manage Dependency panel), instead of hand-editing `setup_Chi_sepnet_environment.m` or a hardcoded path
* New `download_HEIDI_toolbox.m` script to automatically download the HEIDI package and register it in `SpecifyToolboxesDirectory.m`; the package itself is hosted as a GitHub Release asset on the SEPIA repo (tag `heidi-sepiaready-v1`, kept separate from SEPIA's own version tags)
* New `download_toolboxes.m` script to check/download FANSI, HEIDI and Tensor-MPPCA in one go, instead of running each toolbox's own setup script separately
* New `setup_sepia.m` script that auto-creates a machine-local `SpecifyToolboxesDirectory.m` from `SpecifyToolboxesDirectory.template.m` the first time it's missing

**Preprocessing**
* New Tensor-MPPCA denoising option (automatically downloads the required external toolbox on first use)
* New upsampling option for the phase/magnitude data prior to processing
* Added STI-Suite's V-SHARP 2D as a background field removal method for multi-slice/2D EPI acquisitions
* New "no unwrapping" option when only field mapping is required (e.g. for functional QSM)

**Masking**
* New two-pass masking option, and a new, generalised mask refinement pipeline (shared between the I/O panel's "Refine mask" option and the QSM panel's two-pass masking)
* Brain extraction is no longer limited to FSL's BET: a new method dropdown adds Otsu's-method thresholding and (when FreeSurfer's `mri_synthstrip` is available) SynthStrip and SynthStrip (no CSF)
* Built-in V-SHARP: fixed a bug where the k-space deconvolution step was missing, causing incomplete background field removal; kernel radius is now specified in mm instead of voxels (and supports anisotropic voxel sizes)
* Fixed several compatibility issues between two-pass masking and the mask refinement pipeline (BIDS directory input, QSM/mask-refinement wrapper argument handling)

**R2\* handling**
* The R2* map is now computed once per pipeline run and reused across the mask refinement, unreliable-voxel exclusion, and QSM CSF-masking steps (previously recomputed redundantly); it is only recomputed automatically if the data is subsequently denoised or upsampled

**Configuration & GUI**
* SEPIA can now parse `sepia_config*.m` pipeline configuration files and extract the algorithm parameters directly, storing them in the GUI figure handle
* Updated config-file parsing to keep up with newer pipeline configuration files (e.g. the HEIDI-related fields added this release)
* Various GUI bug fixes for loading saved configuration files (e.g. NDI's GPU option, VSHARP/FANSI parameters)
* GUI default method selection is now toolbox-availability-aware for the total field/phase unwrapping, background field removal, and QSM dipole inversion steps, following a consensus-informed priority chain per step (e.g. QSM defaults to FANSI → MEDI → LSQR+HEIDI → TKD, whichever is actually installed); the background field removal "remove residual B1 field" default (3D Polynomial / None) now automatically follows whichever BFR method is selected
* Added a dark theme for the GUI

**BIDS / I/O**
* Added support for reading multiple volumes per echo in BIDS-formatted data (e.g. functional QSM)
* Fixed echo-tag (`_echo-##_`) parsing to work regardless of zero-padding used in the echo number
* Fixed an undefined input-filename-cell error when auto-detecting a BIDS directory containing single-volume-per-echo data
* Pipeline outputs now include BIDS-Derivatives-style JSON sidecars (units, source files, algorithm parameters) alongside the NIfTI files, plus a `dataset_description.json` at the output root
* Output NIfTI extension (`.nii` vs `.nii.gz`) is now detected from the input data instead of always being forced to `.nii.gz`
* Fixed output filenames ending up with two `desc-` BIDS entities when the output prefix already contained one (e.g. from a previous processing stage); it is now merged with SEPIA's own output-type label instead, chained in actual processing order (e.g. denoised → upsampled)
* Renamed the paramagnetic/diamagnetic susceptibility map outputs from the non-standard `ChiParamap`/`ChiDiamap` suffixes to the BIDS-valid `desc-paramagnetic_Chimap`/`desc-diamagnetic_Chimap`
* Fixed `sepiaIO` not resolving relative input/mask/output paths against the working directory before processing (could break since SEPIA changes its current directory internally mid-run)
* `sepiaIO` now always generates a fresh `sepia_config` file for a run, instead of reusing one already present in the output directory

**Segmentation & analysis**
* Automatic contrast matching, quick nonlinear registration using a dilated subcortical mask, label-based registration, and CSV statistics export added to the MuSus-100/CIT168/AHEAD atlas-based segmentation tools
* Chimap values can now be exported to a CSF file after segmentation

**Bug fixes**
* Fixed `get_set_qsm_ndi.m` erroring when loading a saved configuration file
* Fixed a bug in R2* NLLS mapping
* Fixed direct file loads (e.g. user-supplied R2*/R2 maps in the Chi-separation wrapper) bypassing the odd-matrix-size zero-padding step - they now go through the same loading path as other auxiliary data
* Fixed `get_set_qsm_Chi_separation.m` leaving literal quotation marks in the R2*/R2 edit fields when reloading a saved configuration file, and not re-triggering the solver dropdown's callback on load (so field enable/disable state and the Dr default could be left mismatched with the loaded solver)
* Fixed `check_and_set_SEPIA_header_data.m` dropping the `r2s` field (and its `availableFileList` entry) when passing header/extra data through

**Housekeeping**
* `SpecifyToolboxesDirectory.m` and `SpecifyAtlasDirectory.m` are no longer tracked in git (now machine-specific and gitignored; see their `.template.m` files and the new `setup_sepia.m`)
* Removed a large set of unused/deprecated legacy wrapper and parser files (e.g. the deprecated `parse_varargin_*` argument parsers, the GPU-prototype `cuBackgroundRemovalMacro.m`/`cuQSMMacro.m` wrappers, and a deprecated GUI callback)
* Renamed/reorganised a few internal analysis and R2* utility functions to avoid name clashes with other repositories

### 1.2.2.6 (commit 1790ac6)
* Support read Input/Output information from sepia_config.m 
* Phase DICOM values are rescaled using the max/min values in the data instead of rescale slope/intercept of the NIFTI

### 1.2.2.5 (commit 8630efe)
* Fix the mismatch between SEPIA defined B0 direction and LPCNN when it is not along the z-direction
* Fix the shared library issue when using ROMEO with latest versions of Matlab on Linux (see [here](https://github.com/korbinian90/ROMEO))
* Allow user to define atlases' directory paths

### 1.2.2.4 (commit 9083249)
* Fix bug when importing SEPIA pipeline configuration files (sepia_config.m) to the GUI for using VSHARP and FANSI

### 1.2.2.3 (commit efde35b)
* Fix bug when using BIDS compatible directory input where magnitude images did not utilise the rescale slope and intercept to obtain the true values for R2* mapping

### 1.2.2.2 (commit e53fd99)
* Fix bug when using BIDS compatible directory input where magnitude images did not utilise the rescale slope and intercept to obtain the true values for QSM

### 1.2.2.1 (commit 1f04298)
* Fix bug when using optimum weight total field computation with odd matrix size data

### 1.2.2 (commit d6bb60e)
* Fix bug for non-double type input for MATLAB's strel function
* Make sure all holes inside the ROI mask are filled after the background field removal step
* ROI (brain) mask is applied on the fieldmap regardless of what method is chosen

### 1.2.1.1 (commit 941cd5b)
* Enable option of GPU processing for FANSI and NDI

### 1.2.1 (commit 190dd44)
* Fix bug for data with odd-number matrix size
* Fix bug for missing file when using R2* mapping with NLLS algorithm

### 1.2 (commit d2f54a3)
* Support several deep learning based methods (BFRnet, xQSM, QSMnet+ and LP-CNN) on Linux
* Support atlas-based subcortical structure segmentation (CIT168 Reinforcement learning atlas, MuSus-100 and AHEAD) on Linux and Mac
* Integrate R2* mapping toolbox into SEPIA
* New function to further refine brain mask by thresholding high R2* voxels on brain edges
* When magnitude image is used for NDI, the image will be normalised by the intensity of the 99th percentile of the masked voxels instead of the maximum to improve robustness

Please visit the documentation website for more info regarding the newly supported methods and functions.

### 1.1.1 (commit a7680bb)
* ROMEO is now packaged together with CLEAR-SWI. To accompany these changes, ROMEO_HOME is renamed to MRITOOLS_HOME
* Supported CLEAR-SWI
* Fixed bug: bipolar readout correction implementation in full processing pipeline is different from the one in Phase unwrapping standalone 
* Added GPU compatibility of NDI
* Fixed bug for NDI (M^2 is now used instead of M as weights)
* Added functionality to remove brain mask edge **before** backfround field removal step.

### 1.1.0 (commit 9ffe0e2)
* New backend architecture for SWI/SMWI algorithms which supports add-on feature like QSM processing 
* Better compatibility with ROMEO
* New implementation of bipolar readout phase offset correction (from which no phase unwrapping is required)
* Provide bipolar readout phase offset estimation as an output
* New implementation on incorporating mono-exponential fitting residual to weighting map generation
* Experimental support to export GE real|imaginary image to phase image

### 1.0.1 (commit 3a2b387)
* Fixed bug when phase NIfTI is in wrapped range with non-unity rescale slope (e.g. from Philips' scanners)
* Updated function performing phase conversion from arbitary DICOM values to radian (could result in minor numerical differences compared to previous versions if the input phase NIfTI not in radian)
* Several other minor bugs fixed

### 1.0.0 (commit 8e35aee)
* Support ROMEO as total field computation and phase unwrapping method
* Support MRI susceptibility calculation methods for QSM dipole field inversion
* Support FANSI v3.0 (note that the algorithm parameters are adapted for this version)
* Improve BIDS compartibility with SEPIA
* Update output filenames in accordance with BIDS format 
* Improve the comparability of weighting maps across different datasets and methods

### 0.8.1.1 (commit 52dd20b)
* Fixed bug when using single-echo dataset
* Fixed bug when input phase data in unit of radian with single datatype

### 0.8.1 (commit c78247d)
* Log file and error message file are now paired (last 15 digits in the extension) instead of sorting in simple numerical order
* Log file and error message file are now supported in both GUI and command-based operations (when using ``sepiaIO``)
* When running SEPIA, the current directory will temporaily move to the output directory to avoid overwriting temporary files if multiple processings happen simultaneously
* A SEPIA pipeline configuration file will be automatically generated using ``sepiaIO`` is the output directory does not have any existing configuration file. This would be useful to look up the pipeline used to produce the results when using command-based operationn.
* Bug fix when running FANSI (details [here](https://github.com/kschan0214/sepia/issues/8))
* Bug fix when getting B0 direction from Sagittal or Coronal acquisition (details [here](https://github.com/kschan0214/sepia/issues/10))
* Bug fix when running QSM standalone with magnitude image for regularisation (details [here](https://github.com/kschan0214/sepia/issues/9))
* Bug fix when running MEDI with zeropadding option is not equal to zero
* (For developer) Improved readiility of how the data are loaded in SEPIA, which could make better BIDS compartibility in the future

### 0.8.0 (commmit b4255d8) 
* New layout for input/output panel for data selection
* New output config file, log file and error message file
* New feature to load parameters in config file to the GUI
* New option to save unwrapped echo phase
* New option to exlcude unreliable voxels
* New option to select reference tissue for QSM normalisation
* Support the lastest version of MEDI toolbox (Jan 15, 2020)
* Support bipolar readout correction for total field recovery with MEDI's non-linear fitting algorithm
* Support extra brain extraction (FSL's BET) parameters from MEDI toolbox
* New 'percentage' option for MEDI+0 algorithm
* Support the lastest version of FANSI toolbox (commit dc68c306)
* New option to use [weak harmonic regularisation](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27483) with FANSI
* Support [nonlinear dipole inversion (NDI)](https://github.com/polakd/NDI_Toolbox) as external library
* Support [SEGUE](https://xip.uclb.com/i/software/SEGUE.html) as external library

**Please upload the MEDI toolbox (Jan 15, 2020) and FANSI toolbox (commit dc68c306) to the lastest version for the best performance.**

### 0.7.3 (commmit 68c53bc)

* Support [nonlinear dipole inversion (NDI)](https://github.com/polakd/NDI_Toolbox) as external library
* Support [SEGUE](https://xip.uclb.com/i/software/SEGUE.html) as external library

### 0.7.2 (commmit bf020ce))  
* Support single-echo dataset
* Bug fix with odd-number matrix dimension by zero-padding
* Offload unuse variables to reduce memory usage
* Bug fix for reading NIfTI when the rescale slope and intercept are not 1 and 0

### 0.7.1 (commmit dc51fbe)
* Support simple susceptibility weighted imaging (SWI) and susceptibility map weighted imaging (SMWI) as part of the GUI
* resolved loading/saving NIfTI issue related to 0.7.0 update
* DICOM input is deprecated: the only possible input is NIfTI data
* fixed bug when running MEDI with CSF regularisation
* fixed bug for single echo SWI
* now support automatic magnitude and phase images detection with name containing string "mag" for magnitude image and "ph" for phase image  
* fixed global phase offset with graph-cut phase unwrapping

### 0.7.0 (commmit e66d8e4)
* redesigned log file format; the algorithms and parameters being used are much clearer and neat than before (previous log file cannot work in this version)
* resolved '.nii.nii' issue when using STI suite algorithms
* resolved no. of iterations with FANSI does not change issue
* resolved problematic QSM results with FANSI when an input matrix is an odd number
* resolved excluded unreliable voxels issue when 3D best path algorithm doesn't work
* improved build-in VSHARP results when there are masked voxels on the image edges
* added image erosion function for background field removal algorithms
* get header function is now compatible with the JSON files generated by dcm2niix and dicm2nii

### 0.6.0 (commmit 1c27dc4)  
* updated diretcory structure
* added options to select individual files  
