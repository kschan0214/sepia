# qsm_hub
====================================================================================================
<h1>Introduction</h1>
--------------------------------------------------

Welcome to the beta version of qsm_hub. qsm_hub is a graphical user interface for quantitative
susceptibility mapping (QSM) processing in MATLAB.

It is built based on three toolboxes including [MEDI](http://weill.cornell.edu/mri/pages/qsm.html) (version 2017-11-06), [STI Suite](https://people.eecs.berkeley.edu/~chunlei.liu/software.html) (version 3.0)
and [FANSI](https://www.martinos.org/~berkin/software.html) (version 2.0).

qsm_hub serves with two purposes:
1. acts as a hub to allow user to choose different QSM (pre)processing methods.
2. allows fast method parameter adjustments with the GUI

Instead of focusing providing processing methods of its own (indeed there are some methods
implemented by me), large amount of effort was paid to create wrappers/macros to provide interface
to incorporate different methods.

qsm_hub provides a tool to produce QSM (and related) map easily. It is particularly suitable for
searching the best processing methods for your multi-echo GRE (mGRE) data.

Once you find the best setting of the processing pipeline (and more familiar with qsm_hub), you
might wish to dive into the wrapper functions for batch processing.

Standard QSM data processing usually involves the following procedures:
0. convert DICOM phase values to wrapped phase values
1. Phase unwrapping and total field recovery
2. Background field removal
3. QSM (sometimes people just call it dipole field inversion)

qsm_hub provides 4 standalone for the above procedures:
1. QSMHub (One-stop QSM processing): one-stop platform from loading the mGRE data(either NIfTI or
	DICOM) to generating susceptibility map
2. Phase unwrapping: standalone to convert complex-valued mGRE data (DICOM or NIfTI) to unwrapped
	total field map
3. Background field removal: standalone to remove background field contribution from a total field
	map to produce a local field map
4. QSM: standalone to map magnetic susceptibility source from a local field map

Apparently this toolbox is still in development, so you should expect to encounter some bugs.

If you have any question or you would like to report bug please feel free to contact me
k.chan@donders.ru.nl (Kwok-Shing Chan).

Have fun!

Kwok (2018-05-07)

----------------------------------------------------------------------------------------------------
<h1>Installation</h1>

To use the qsm_hub, just add the directory containing qsm_hub.m into your MATLAB PATH

This can be done by: 'Set Path' -> 'Add Folder' -> *your directory of qsm_hub.m* -> 'Save'

or

with MATLAB's command: addpath('/your/qsm_hub/directory')

----------------------------------------------------------------------------------------------------
<h1>Compatibility</h1>

I try to make qsm_hub compatible with as many MATLAB version as I could. In principle, it works fine
with MATLAB R2016b or earlier. However, most of the methods should also be working with MATLAB
R2017a or later, except 'LBV' of background field removal method.

----------------------------------------------------------------------------------------------------
<h1>Working with NIfTI files and qsm_hub header</h1>


If you prefer working with NIfTI files (as I do), I suggest using MRIConvert
(https://lcni.uoregon.edu/downloads/mriconvert) to convert your mGRE data with 'Option' ->
'Save multivolumes series as 4D files' checked. In this way your mGRE data will be stored as 4D FSL
NIfTI data that needed to be the input of qsm_hub. A text file will also be generated with the NIfTI
data. If your put this text file in your qsm_hub input directory, the actual echo times stated in
this text file will also be read to generate a synthetic qsm_hub header.

qsm_hub requires a special header file in mat format that stores some information for QSM recon. If
your input are DICOMs, this header file will be generated automatically. For NIfTI data, if this
file is absent, qsm_hub will generate a synthetic header based on the information from the NIfTI
files. However, for some information such as 'TE' and 'delta_TE', the NIfTI header usually doesn't
contain them and there some predefined values will be set. The synthetic header may affect the
susceptibility range of the QSM value but not the qualitative accessment QSM (basiclly the QSMs are
the same but in different (arbitary) scale)

Altenatively you can create the qsm_hub header file in your own way, but please make sure that the
header filename contains the string 'header' and the file contains the following variables:
	-	'B0'			: magnetic field strength, in Telsa (e.g. B0=3 % 3T)
	-	'B0_dir'		: main magnetic field direction, [x,y,z] (e.g. B0_dir=[0,0,1])
	-	'CF'			: imaging frequency, in Hz (e.g. CF=3*42.58*1e6 %water 1H at 3T)
	-	'TE' 			: echo times, in s (e.g. TE=[te1,te2,te3,te4,te5])
	-	'delta_TE'		: echo spaing, in s (e.g. delta_TE=TE(2)-TE(1))
	-	'matrixSize'	: image matrix size (e.g. matrixSize=size(img))
	-	'voxelSize'		: spatial resolution of the data (e.g. voxelSize=[2,2,2] % 2 mm isotropic)
I encourage to use SyntheticQSMHubHeader.m to get most of the information from NIfTI header such as
'B0_dir', 'matrixSize' and 'voxelSize', and readTEfromText.m to get echo time information.

----------------------------------------------------------------------------------------------------
<h1>Checking QSM result</h1>


Every now and again the values of QSM might be inverted. To check if it is the case you could look
into the QSM map. The deep gray matter structure should be appeared as bright (positive) and white
matter should be dark (negative value). If it appears in the opposite way then you might invert the
result by multiply the map with -1, i.e. QSM_corr = -QSM.

----------------------------------------------------------------------------------------------------
<h1>Standalone description</h1>


----------------------------------
|QSMHub (One-stop QSM processing)|
----------------------------------
	I/O planel
--------------------------------------------------
- Input:
(Option 1) 	Directory contains all and only mGRE DICOM data, including both magnitude and phase
			images
(Option 2) 	Directory contains both magnitude and phase 4D-NIfTI data([x,y,slice,time]) and qsmhub
			header mat file (e.g. qsmhub_header.mat, optional), the magnitude data filename must
			contain the string 'magn' (e.g. 'magn.nii.gz') and the phase data filename must contain
			the string 'phase' (e.g. 'phase.nii.gz')
- Output:
	Directory to store all qsm_hub output (default: /input/dir/output)
- FSL brain extraction (optional):
	Simple FSL's BET script, which can roughly extract a brain yet not too accurate.
	User is recommended to perform accurate brain extraction outside the qsm_hub app (e.g. using
	FSL's bet which allows fine tuning for more accurate result)
- Brain mask file:
	NIfTI file that brain voxel > 0.
	QSM relies on brain mask thorough all processing steps.
	If the input directory also contains NIfTI file that contains string 'mask' in filename (e.g.
	'mask.nii.gz'), it will be read automatically so you might skip this steps

[Example]
A standard input directory contains the following files:
- 	magn.nii.gz			(4D real image)
- 	phase.nii.gz		(4D real image)
- 	mask.nii.gz 		(3D mask image)
- 	qsmhub_header.mat	(.mat file)

	Total field recovery and phase unwrapping
--------------------------------------------------
- Method:
	(1) Laplacian: 				very reliable method for phase unwrapping yet the output values
								are not accurate
	(2) Laplacian STI suite: 	Laplacian unwrapping implementation from STI Suite v3.0
	(3) Jena: 					very robust region growing method yet only works in the DCCN cluster
								(recommended in the cluster)
	(4) Region growing: 		MEDI toolbox implementation, might not work well with DICOM phase
								data (using offline recon data works pretty well)
	(5) Graphcut: 				graph-cut algorithm (not recommended)
- Bipolar readout eddy current correction:
	enable to correct the phase inconsistence between odd and even echoes and gradient field
	contribution by eddy current effect due to bipolar readout
- Exclude unreliable voxels, Threshold:
	enable to exclude low SNR voxels that can create strong artefact in susceptibility map (you may
	check with 'squirrel_fieldmapSD.nii.gz' to adjust the threshold)
[Output]
- 	squirrel_phase_EC.nii.gz 		(eddy current corrected phase, optional)
- 	squirrel_magn_EC.nii.gz			(eddy current corrected magnitude, optional)
- 	squirrel_totalField.nii.gz 		(unwrapped total (background+local) field, in Hz)
- 	squirrel_fieldMapSD.nii.gz 		(normalised field map standaed deviation)

	Background field removal
--------------------------------------------------
- Method:
	(1)	LBV: 				Laplacian boundary value approach to removal background field
							(recommended)
	(2) PDF: 				Project onto dipole field
	(3) RESHARP: 			Regularised SHARP
	(4) SHARP: 				Sophisticated harmonic artefact reduction for phase data
	(5) VSHARP STI suite: 	STI suite v3.0 variable-kernel SHARP (recommended)
	(6) VSHARP
	(7) iHARPERELLA: 		(not recommended)
- Refine local field by 4th order 3D polynomial fit
	enable to remove residual B1(+ & -) contribution in the local field
[Output]
- 	squirrel_localField.nii.gz 		(background-field-removed local (or tissue) field, in Hz)
- 	squirrel_mask_final.nii.gz		(final brain mask, might be eroded depended on background field
									removal algorithms and 'exclude unrelaiable voxels' threshold value)

	QSM
--------------------------------------------------
- Method:
	(1) TKD: 					Thresholded k-space division
	(2) Closed-form solution: 	closed-form solution with L2-norm regularisation
	(3) STI suite iLSQR: 		STI suite v3.0 implementation of iterative LSQR approach
	(4) iLSQR
	(5) FANSI: 					Fast algorithm for nonlinear susceptibility inversion (recommended)
	(6) Star: 					STI suite v3.0 Star-QSM (recommended)
	(7) MEDI: 					Morphology enabled dipole inversion (MEDI+0) (pretty good but slow)
[Output]
- 	squirrel_QSM.nii.gz 			(quantitative susceptibility map, in ppm)

------------------
|Phase unwrapping|
------------------
	I/O planel
--------------------------------------------------
- Input:
(Option 1) 	Directory contains all and only mGRE DICOM data, including both magnitude and phase
			images
(Option 2) 	Directory contains both magnitude and phase 4D-NIfTI data([x,y,slice,time]) and qsmhub
			header mat file (e.g. qsmhub_header.mat, optional), the magnitude data filename must
			contain the string 'magn' (e.g. 'magn.nii.gz') and the phase data filename must contain
			the string 'phase' (e.g. 'phase.nii.gz')
- Output:
	Directory to store all qsm_hub output (default: /input/dir/output)
- FSL brain extraction (optional):
	Simple FSL's BET script, which can roughly extract a brain yet not too accurate.
	User is recommended to perform accurate brain extraction outside the qsm_hub app (e.g. using
	FSL's bet which allows fine tuning for more accurate result)
- Brain mask file:
	NIfTI file that brain voxel > 0.
	QSM relies on brain mask thorough all processing steps.
	If the input directory also contains NIfTI file that contains string 'mask' in filename (e.g.
	'mask.nii.gz'), it will be read automatically so you might skip this steps

	Total field recovery and phase unwrapping
--------------------------------------------------
- Method:
	(1) Laplacian: 				very reliable method for phase unwrapping yet the output values
								are not accurate
	(2) Laplacian STI suite: 	Laplacian unwrapping implementation from STI Suite v3.0
	(3) Jena: 					very robust region growing method yet only works in the DCCN cluster
								(recommended in the cluster)
	(4) Region growing: 		MEDI toolbox implementation, might not work well with DICOM phase
								data (using offline recon data works pretty well)
	(5) Graphcut: 				graph-cut algorithm (not recommended)
- Bipolar readout eddy current correction:
	enable to correct the phase inconsistence between odd and even echoes and gradient field
	contribution by eddy current effect due to bipolar readout
- Exclude unreliable voxels, Threshold:
	enable to exclude low SNR voxels that can create strong artefact in susceptibility map (you may
	check with 'squirrel_fieldmapSD.nii.gz' to adjust the threshold)
[Output]
- 	squirrel_phase_EC.nii.gz 		(eddy current corrected phase, optional)
- 	squirrel_magn_EC.nii.gz			(eddy current corrected magnitude, optional)
- 	squirrel_totalField.nii.gz 		(unwrapped total (background+local) field, in Hz)
- 	squirrel_fieldMapSD.nii.gz 		(normalised field map standaed deviation)

--------------------------
|Background field removal|
--------------------------
	I/O planel
--------------------------------------------------
- Input:
	Directory contains total field map, (optional) field map standard deviation map NIfTI data and qsmhub
	header mat file (e.g. qsmhub_header.mat, optional).
	The total field map data filename must contain the string 'totalfield' (e.g.
	'squirrel_totalField.nii.gz') and the field map standard deviation data filename must contain the
	string 'fieldmapsd' (e.g. 'squirrel_fieldMapSD.nii.gz').
- Output:
	Directory to store local field map and final brain mask (default: /input/dir/output)
- Brain mask file:
	NIfTI file that brain voxel > 0. QSM relies on brain mask thorough all processing steps. If the input
	directory also contains NIfTI file that contains string 'mask' in filename (e.g. 'mask.nii.gz'), it
	doesn't need to be specified and will be read automatically
[Example]
A standard input directory contains the following files:
- 	totalField.nii.gz	(unwrapped total field map)
- 	fieldMapSD.nii.gz	(normalised field map standard deviation, optional)
- 	mask.nii.gz 		(3D mask image)
- 	qsmhub_header.mat 	(.mat file)

	Background field removal
--------------------------------------------------
- Method:
	(1) LBV: Laplacian boundary value approach to removal background field (recommended)
	(2) PDF: Project onto dipole field
	(3) RESHARP: Regularised SHARP
	(4) SHARP: Sophisticated harmonic artefact reduction for phase data
	(5) VSHARP STI suite: STI suite v3.0 variable-kernel SHARP (recommended)
	(6) VSHARP
	(7) iHARPERELLA: (not recommended)
- Refine local field by 4th order 3D polynomial fit
	enable to remove residual B1(+ & -) contribution in the local field
[Output]
- 	squirrel_localField.nii.gz 	(background-field-removed local (or tissue) field, in Hz)
- 	squirrel_mask_final.nii.gz	(final brain mask, might be eroded depended on background field removal
								algorithms and 'exclude unrelaiable voxels' threshold value)

-----
|QSM|
-----
	I/O planel
--------------------------------------------------
- Input:
	Directory contain local field map and (optional) magnitude mGRE NIfTI data and qsmhub header mat file
	(e.g. qsmhub_header.mat, optional)
	The total field map data filename must contain the string 'totalfield' (e.g. 'squirrel_totalField.nii.gz')
	and the field map standard deviation data filename must contain the string 'fieldmapsd' (e.g.
	'squirrel_fieldMapSD.nii.gz')
- Output:
	Directory to store local field map and final brain mask (default: /input/dir/output)
- Brain mask file:
	NIfTI file that brain voxel > 0. QSM relies on brain mask thorough all processing steps. If the input
	directory also contains NIfTI file that contains string 'mask' in filename (e.g. 'mask.nii.gz'), it doesn't
	need to be specified and will be read automatically
[Example]
A standard input directory contains the following files:
- 	localField.nii.gz 	(3D local field map, in Hz)
- 	magn.nii.gz			(4D real image used to cerate weighting map for some QSM algorithm, [x,y,slice,time],
						optional)
- 	mask_final.nii.gz 	(3D mask image)
- 	qsmhub_header.mat 	(.mat file)

	QSM
--------------------------------------------------
- Method:
	(1) TKD: Thresholded k-space division
	(2) Closed-form solution: closed-form solution with L2-norm regularisation
	(3) STI suite iLSQR: STI suite v3.0 implementation of iterative LSQR approach
	(4) iLSQR
	(5) FANSI: Fast algorithm for nonlinear susceptibility inversion (recommended)
	(6) Star: STI suite v3.0 Star-QSM (recommended)
	(7) MEDI: morphology enabled dipole inversion (MEDI+0) (pretty good but slow)
[Output]
- 	squirrel_QSM.nii.gz 			(quantitative susceptibility map, in ppm)
