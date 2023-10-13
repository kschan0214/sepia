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
- [mritools (ROMEO/CLEARSWI) (v3.5.5)](https://github.com/korbinian90/CompileMRI.jl/releases) (2022-Oct-11: v3.5.6 also passed),
- [MRI Susceptibility Calculation Methods, accessed 12 September 2019](https://xip.uclb.com/product/mri_qsm_tkd).

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

### 1.2.2.5 (current master)
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
