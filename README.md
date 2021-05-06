# SEPIA (SuscEptibility mapping PIpeline tool for phAse images)

![sepia logo](https://sepia-documentation.readthedocs.io/en/latest/_static/logo.png)

## Introduction  

**SEPIA** is a tool providing a graphical user interface to build data processing pipeline of quantitative susceptibility mapping (QSM) in Matlab.

This GUI is built to access the following toolboxes:
[MEDI](http://weill.cornell.edu/mri/pages/qsm.html), 
[STI Suite](https://people.eecs.berkeley.edu/~chunlei.liu/software.html),
[FANSI](https://gitlab.com/cmilovic/FANSI-toolbox),  
[SEGUE](https://xip.uclb.com/i/software/SEGUE.html), and 
[nonlinear dipole inversion (NDI)](https://github.com/polakd/NDI_Toolbox).

SEPIA provides two key features for QSM processing:  
1. mix-and-match methods from different toolboxes to build your own QSM processing pipeline,
2. graphical user interface to easily adjust parameters of different algorithms.

SEPIA is designed to provide a platform for easy access to different QSM processing methods in the field. To achieve this, most of the codes were written for data flow and algorithm parameter control. Through SEPIA, I hope researchers who are not expert in QSM will also be able to use QSM for their research.

**For better readability, the documentation of SEPIA has moved to https://sepia-documentation.readthedocs.io/.**  

## Terms of use
All the original codes and methods developed for **SEPIA** are under MIT license. You can check [the license file](https://github.com/kschan0214/Sepia/blob/master/LICENSE) for more information. For the terms of use of the toolboxes related to this work, their own license applied and please check the corresponding license file(s) in each toolbox for more information. 

If you use SEPIA in your research, please cite the following article:

[Chan, K.-S., Marques, J.P., 2021. SEPIA—Susceptibility mapping pipeline tool for phase images. Neuroimage 227, 117611.](https://doi.org/10.1016/j.neuroimage.2020.117611)  

As well as any related papers in your processing pipeline. 

If you have any question or you would like to provide suggestion to improve this toolbox/report bug(s) please raise an issue in the [github page](https://github.com/kschan0214/sepia/issues).


## Update notes  

For full update log, please visit https://sepia-documentation.readthedocs.io/.

### Future release
* Support FANSI v2.0
* Support [iterative Tikhonov regularisation for QSM dipole inversion](https://xip.uclb.com/i/software/mri_qsm_tkd.html) 
* Better compartibility with BIDS format data

### 0.8.1.1 (current master)
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
