# qsm_hub

Welcome to the beta version of qsm_hub. qsm_hub is a graphical user interface for quantitative susceptibility mapping (QSM) processing. It acts as a hub to allow user to choose most of the QSM (pre)processing methods. Standard QSM data processing steps involve the following procedures:
[(0) convert DICOM phase values to wrapped phase values]
(1) Phase unwrapping and total field recovery
(2) Background field removal
(3) QSM (a.k.a. dipole field inversion)

qsm_hub provides 4 standalone for the above procedures:
(1) QSMHub (One-stop QSM processing): one-stop platform from loading the (m)GRE data(either NIfTI or DICOM) to generating susceptibility map
(2) Phase unwrapping: standalone to convert complex-valued (m)GRE data (DICOM or NIfTI) to unwrapped total field map
(3) Background field removal: standalone to remove background field contribution from a total field map to produce a local field map
(4) QSM: standalone to map magnetic susceptibility source from a local field map

QSMHub:
I/O option -
Input: 
(Option 1) Directory contains all and only mGRE DICOM data, including both magnitude and phase images
(Option 2) Directory contain both magnitude and phase 4D-NIfTI data([x,y,slice,time]) and qsmhub header mat file (e.g. qsmhub_header.mat, optional), the magnitude data filename must contain the string 'magn' (e.g. 'magn.nii.gz') and the phase data filename must contain the string 'phase' (e.g. 'phase.nii.gz')
Output:
Directory to store all qsm_hub output (default: /input/dir/output)
FSL brain extraction (optional):
Simple FSL's BET script, which can roughly extract a brain yet not too accurate. User is recommended to perform accurate brain extraction outside qsm_hub app (e.g. using FSL's bet which allows fine tune with hyper parameters fr better accuracy)
Brain mask file:
NIfTI file that brain voxel > 0. QSM relies on brain mask thorough all processing steps. If the input directory also contains NIfTI file that contains string 'mask' in filename (e.g. 'mask.nii.gz'), it doesn't need to be specified and will be read automatically

A standard input directory contains the following files:
- magn.nii.gz	(4D real image)
- phase.nii.gz	(4D real image)
- mask.nii.gz 	(3D mask image)
- qsmhub_header.mat (.mat file contains the following variables 'B0' (field strength), 'B0_dir' (main magnetic field direction), 'CF' (imaging frequency), 'TE' (all echo times), 'delta_TE' (echo spacing, a.k.a. echo time), 'matrixSize' (matrix size of 3D image), 'voxelSize' (spatial resolution). This is a file that only being needed for NIfTI files input: for DICOM data it will be automatically generated. If this file is absent, qsm_hub will generate a synthetic header based on the information from the NIfTI files. However, for some information such as 'TE' and 'delta_TE', the NIfTI header usually doesn't contain these information. The synthetic header may affect the susceptibility range of the QSM value)

Total field recovery and phase unwrapping - 
Method:
(1) Laplacian: very reliable method for phase unwrapping yet the output values are not accurate
(2) Laplacian STI suite: Laplacian unwrapping implementation from STI Suite v3.0
(3) Jena: very robust region growing method yet only works in the cluster (recommended in the cluster)
(4) Region growing: MEDI toolbox implementation, might not work well with DICOM phase data (using offline recon data works pretty well)
(5) Graphcut: graph-cut algorithm (not recommended)
Bipolar readout eddy current correction:
enable to correct the phase inconsistence between odd and even echoes and gradient field contribution by eddy current effect due to bipolar readout
Exclude unreliable voxels, Threshold:
enable to exclude low SNR voxels that can create strong artefact in susceptibility map(you may check with 'squirrel_fieldmapSD.nii.gz' to adjust the threshold)

Background field removal - 
Method:
(1) LBV: Laplacian boundary value approach to removal background field (recommended)
(2) PDF: Project onto dipole field
(3) RESHARP: Regularised SHARP
(4) SHARP: Sophisticated harmonic artefact reduction for phase data 
(5) VSHARP STI suite: STI suite v3.0 variable-kernel SHARP (recommended)
(6) VSHARP
(7) iHARPERELLA: (not recommended)
Refine local field by 4th order 3D polynomial fit
enable to remove residual B1(+ & -) contribution in the local field

QSM - 
Method:
(1) TKD: Thresholded k-space division
(2) Closed-form solution: closed-form solution with L2-norm regularisation 
(3) STI suite iLSQR: STI suite v3.0 implementation of iterative LSQR approach
(4) iLSQR
(5) FANSI: Fast algorithm for nonlinear susceptibility inversion (recommended)
(6) Star: STI suite v3.0 Star-QSM (recommended)
(7) MEDI: morphology enabled dipole inversion (MEDI+0) (pretty good but slow)
  
If you have any question or you would like to report bug please feel free to contact me k.chan@donders.ru.nl (Kwok-Shing Chan).

Have fun!

Kwok