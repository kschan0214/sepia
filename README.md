# QSM Hub (work in progress)
==================================================

## Introduction
--------------------------------------------------

Welcome to the beta version of **QSM Hub**. **QSM Hub** is a tool providing a graphical user 
interface to build data processing pipeline of quantitative susceptibility mapping (QSM) in Matlab.

This GUI is built based on three toolboxes including [MEDI](http://weill.cornell.edu/mri/pages/qsm.html) 
(version 2017-11-06), [STI Suite](https://people.eecs.berkeley.edu/~chunlei.liu/software.html) (version 3.0)
and [FANSI](https://www.martinos.org/~berkin/software.html) (version 2.0).

**QSM Hub** provides two key features for QSM processing:

1. mix-and-match methods from different toolboxes to build your own QSM processing pipeline
2. graphical user interface to easily adjust parameters of different algorithms

The objective of **QSM Hub** is to provide a platform for easy access of different QSM processing 
method in the field. To achieve this, most of the codes were written for data flow control and 
algorithm parameter control. Through **QSM Hub** I hope researchers which are not coming from the field 
would also be able to use QSM for their research. For those more experienced researchers familiar with 
QSM, I also hope this tool can simplify your workflow. 

A Standard QSM data processing piepline usually consists of the following procedures:

0. convert DICOM phase values (4096,4095) to wrapped phase values (-pi,pi) (may not be neccessary)  
1. phase unwrapping and total field recovery  
2. background field removal  
3. quantitative susceptibility mapping  

**QSM Hub** provides 4 standalone applications for the above procedures:  

1. **QSMHub (One-stop QSM processing)**  
	one-stop processing platform from loading the GRE data(either NIfTI or DICOM) to generate a susceptibility map  
2. **Phase unwrapping**  
	standalone that converts complex-valued mGRE data (DICOM or NIfTI) to unwrapped total field map  
3. **Background field removal**  
	standalone that removes background field contribution from a total field map to produce a local field map  
4. **QSM**  
	standalone that maps magnetic susceptibility sources from a local field map  

This toolbox is still under development, so you may encounter some bugs.

**QSM Hub** is freely available for educational and research purposes. Commerical or clinical use is not 
permitted. When you use **QSM Hub** in your research, please cite the related papers of the methods 
in your processing pipeline. Details of each processing method is provided with the hyperlinks in this
document. Clicking the hyperlinks of the methods below will redirect you to the 
paper webiste or you can also find the full details in the Referencing section below. 

If you have any question or you would like to provide suggestion to improve this toolbox/report 
bug(s) please feel free to contact me k.chan@donders.ru.nl (Kwok-Shing Chan).

Kwok  
2018-05-30

----------------------------------------------------------------------------------------------------

## Installation
--------------------------------------------------

To use **QSM Hub**, just add the directory containing qsm_hub.m into your MATLAB PATH

This can be done by:  
'Set Path' -> 'Add Folder' -> /your/qsm_hub/directory/ -> 'Save'  
  
or
  
with MATLAB's command: `addpath('/your/qsm_hub/directory')`  

Then, you can start the GUI by entering `qsm_hub` in the MATLAB's command window.  

----------------------------------------------------------------------------------------------------

## Compatibility
--------------------------------------------------

In principle, **QSM Hub** works fine with MATLAB R2016b or earlier. However, most of the methods 
should also be working with MATLAB R2017a or later, except 'LBV' of background field removal method.  

----------------------------------------------------------------------------------------------------

## Input data format and qsm_hub header
--------------------------------------------------

**QSM Hub** provides two main options for data input:  

1. DICOM images    
2. NIfTI images (.nii and .nii.gz)  

You can find the specific requirements of the input files in the 'standlone description' section.

If you prefer starting with the unprocessed data in NIfTI format , **QSM Hub** expects the input to 
be 4D (row,col,slice,echo). I suggest using 
[MRIConvert](https://lcni.uoregon.edu/downloads/mriconvert) to convert your mGRE data with the 
following setting checked:  
'Option' -> 'Save multivolumes series as 4D files'  
In this way your mGRE data will be stored as 4D FSL NIfTI data that is a valid input of **QSM Hub**. 

**QSM Hub** requires a specific MAT-file (.mat) that stores some header information of the input 
data. The header file contains the following variables (case-sensitive):
````
	'B0'			: magnetic field strength, in Tesla (e.g. B0=3; % 3T)
	'B0_dir'		: main magnetic field direction, [x,y,z] (e.g. B0_dir=[0,0,1];)
	'CF'			: imaging frequency, in Hz (e.g. CF=3*42.58*1e6; %water 1H at 3T)
	'TE' 			: echo times, in s (e.g. TE=[te1,te2,te3,te4,te5];)
	'delta_TE'		: echo spacing, in s (e.g. delta_TE=TE(2)-TE(1);)
	'matrixSize'	: image matrix size (e.g. matrixSize=size(img);)
	'voxelSize'		: spatial resolution of the data (e.g. voxelSize=[2,2,2]; % 2 mm isotropic)
````  
If your input is DICOM images, this file will be generated automatically with **QSM Hub**. If your 
input is NIfTI images, you can use the utility function provided called 'Get qsm_hub header' in the 
'Utility' tab of the GUI to generate this header file. Once this file is generated, make sure its 
name contains the string 'header' (e.g. 'qsmhub_header.mat') and put it in the same directory of 
your **QSM Hub** input data. 

**Tip**:
If you use MRIConvert to convert DICOM images to NIfTI format, a text file will also be generated 
alongside with the NIfTI data. This text file contains some basic header information including the 
echo times. When you use the 'Get qsm_hub header' utility function, you can load this text file as 
the 'TE file' so that the correct echo time information can be extracted and stored in the MAT-file 
header.

----------------------------------------------------------------------------------------------------

## Checking QSM result
--------------------------------------------------

Sometimes the values of QSM might be inverted. This may be due to the way of how the phase 
data being read. To check if this is the case you could look into the QSM map. The deep gray matter 
structure should be appeared as bright (positive value) 
If it appears in the opposite way then you might invert the result by multiply the map with -1, i.e.  

`QSM_corr = -QSM; `

or checked the 'Invert phase' option in the GUI.

----------------------------------------------------------------------------------------------------

## GPU computation
--------------------------------------------------

Some QSM processing algorithms utilise iterative approaches (e.g. `pcg`, `lsqr`, etc.) to obatain a 
stable resulting map. To ensure convergence, one can be done is to increase the number of 
iterations but the computation time will also be increasing. **QSM Hub** supports GPU computation in 
Matlab. The GUI function will automatically detect the number of GPU devices available in local 
computer (using `gpuDeviceCount`). If the GPU device is detectable in Matlab, the GUI will provide
an option for user to choose using the GPU feature or not. Based on some preliminary testing, the 
following methods will be benefical (and/or workable) with GPU computation:  

1. Projection onto dipole field (background field removal)  
2. Closed-form L2-regularisation (QSM)  
3. iLSQR (QSM)  
4. FANSI (QSM)  
5. Star (QSM)  
6. MEDI (QSM)  

**Caution**  
The performance of GPU computation is affected by not only the number of GPUs but also the RAM of
the GPU device. In principle, the input data of the methods will only be loaded onto the GPUs when 
it is being called and the GPU device will be reset after the computation. However, the GPU 
implementations could still be broken if the size of input data is large (e.g. large matrix size 
data). The rule of thumb is that if you aim at a small iteration number (e.g. 30-50, except the 
case of using MEDI), the performance of CPU and GPU computation will not have significant 
differences (sometimes GPU may be even slower due to data distribution prior computation).
Nevertheless, if you are trying more iterations (i.e. >100), GPU computation usually gives you some
descent speed improvement.  

This feature is only tested with the following system configuration:  
- macOS High Sierra 10.13.4  
- NVIDIA GeForce GT 750M  
- Matlab R2018a  
If it doesn't work on your system, you can simply disable (unchecked) the GPU option.

----------------------------------------------------------------------------------------------------

## Acknowledgements
--------------------------------------------------

The methods and algorithms provided with **QSM Hub** are contributed by:  

[**MEDI Toolbox**](http://weill.cornell.edu/mri/pages/qsm.html)  
Weill Cornell MRI lab   
- DICOM reader  
- Matlab interface of using FSL's brain extraction  
- Combination of echo phase using non-linear fitting  
- Laplacian phase unwrapping  
- Region growing phase unwrapping  
- Graphcut phase unwrapping  
- Laplacian boundary value (LBV) background field removal  
- Projection onto dipole field (PDF) background field removal  
- Morphology enabled dipole inversion QSM  
- Lateral ventricle masking  

[**STI Suite**](https://people.eecs.berkeley.edu/~chunlei.liu/software.html)  
Hongjiang Wei, PhD - University of California Berkeley, CA, USA  
Wei Li, PhD, - University Texas Health Science Center at San Antonio Research Imaging Institute, San 
Antonio, TX, USA  
Chunlei Liu, PhD - University of California Berkeley, CA, USA    
- Laplacian phase unwrapping  
- VSHARP background field removal  
- iHARPERELLA background field removal  
- iLSQR QSM    
- Star-QSM   

[**FANSI Toolbox**](https://martinos.org/~berkin/software.html)  
Carlos Milovix, PhD - Biomedical Imaging Center at Pontificia Universidad Católica de Chile   
Berkin Bilgic, PhD - Martinos Center for Biomedical Imaging, Harvard Medical School, MA, USA    
Bo Zhao, PhD - Martinos Center for Biomedical Imaging, Harvard Medical School, MA, USA   
Julio Acosta-Cabronero, PhD - Wellcome Trust Centre for Neuroimaging, Institute of Neurology, 
University College London, London, UK, and German Center for Neurodegenerative Diseases (DZNE), 
Magdeburg, Germany   
Cristian Tejos, PhD - Department of Electrical Engineering, Pontificia Universidad Catolica de 
Chile, Santiago, Chile and the Biomedical Imaging Center at Pontificia Universidad Católica de Chile  
- FANSI QSM (including both linear and non-linear solvers, and TV and TGV regularisation options)   

----------------------------------------------------------------------------------------------------

## Referencing
--------------------------------------------------

When you can **QSM Hub** in your research, please cite the method(s) that you used:

### Brain extraction
==================================================  
**FSL bet**  
[Smith, S. M. Fast robust automated brain extraction. Hum. Brain Mapp. 17, 143–155 (2002).
](https://doi.org/10.1002/hbm.10062)

### Phase unwrapping  
==================================================

**Laplacian-based method**   
[Schofield, M. A. & Zhu, Y. Fast phase unwrapping algorithm for interferometric applications. Opt 
Lett 28, 1194–1196 (2003).](https://doi.org/10.1364/OL.28.001194)  

[Li, W., Wu, B. & Liu, C. Quantitative susceptibility mapping of human brain reflects spatial 
variation in tissue composition. 
Neuroimage 55, 1645–1656 (2011).](https://doi.org/10.1016/j.neuroimage.2010.11.088)  

**3D best path**  
[Abdul-Rahman, H. S. et al. Fast and robust three-dimensional best path phase unwrapping algorithm. 
Applied Optics 46, 6623–6635 (2007).](https://doi.org/10.1364/AO.46.006623)  

**Graphcut**    
[Dong, J. et al. Simultaneous phase unwrapping and removal of chemical shift (SPURS) using graph 
cuts: application in quantitative susceptibility mapping. IEEE Transactions on Medical Imaging 34, 
531–540 (2015).](https://doi.org/10.1109/TMI.2014.2361764)    

**Echo phase combination** - **Optimum weights**   
[Robinson, S. D. et al. An illustrated comparison of processing methods for MR phase imaging and QSM: 
combining array coil signals and phase unwrapping. 
NMR Biomed 30, e3601 (2017).](https://doi.org/10.1002/nbm.3601)   

**Echo phase combination** - **MEDI nonlinear fit**  
[Liu, T. et al. Nonlinear formulation of the magnetic field to source relationship for robust 
quantitative susceptibility mapping. 
Magn Reson Med 69, 467–476 (2012).](https://doi.org/10.1002/mrm.24272)  


### Background field removal  
==================================================

**LBV**    
[Zhou, D., Liu, T., Spincemaille, P. & Wang, Y. Background field removal by solving the Laplacian 
boundary value problem. NMR Biomed 27, 312–319 (2014).](https://doi.org/10.1002/nbm.3064)   

**PDF**  
[Liu, T. et al. A novel background field removal method for MRI using projection onto dipole 
fields (PDF). NMR Biomed 24, 1129–1136 (2011).](https://doi.org/10.1002/nbm.1670)    

**RESHARP**    
[Sun, H. & Wilman, A. H. Background field removal using spherical mean value filtering and Tikhonov 
regularization. Magn Reson Med 71, 1151–1157 (2014).](https://doi.org/10.1002/mrm.24765)    

**SHARP**  
[Schweser, F., Deistung, A., Lehr, B. W. & Reichenbach, J. R. Quantitative imaging of intrinsic 
magnetic tissue properties using MRI signal phase: an approach to in vivo brain iron metabolism? 
Neuroimage 54, 2789–2807 (2011).](https://doi.org/10.1016/j.neuroimage.2010.10.070)    

**VSHARP**   
[Li, W., Wu, B. & Liu, C. Quantitative susceptibility mapping of human brain reflects spatial 
variation in tissue composition. 
Neuroimage 55, 1645–1656 (2011).](https://doi.org/10.1016/j.neuroimage.2010.11.088)  

**iHARPERELLA**  
[Li, W., Avram, A. V., Wu, B., Xiao, X. & Liu, C. Integrated Laplacian-based phase unwrapping and 
background phase removal for quantitative susceptibility mapping. 
NMR Biomed 27, 219–227 (2014).](https://doi.org/10.1002/nbm.3056)  

### QSM
==================================================

**TKD**  
[Wharton, S., Schäfer, A. & Bowtell, R. Susceptibility mapping in the human brain using 
threshold-based k-space division. 
Magn Reson Med 63, 1292–1304 (2010).](https://doi.org/10.1002/mrm.22334)  

[Shmueli, K. et al. Magnetic susceptibility mapping of brain tissue in vivo using MRI phase data. 
Magn Reson Med 62, 1510–1522 (2009).](https://doi.org/10.1002/mrm.22135)  

**Closed-form solution**  
[Bilgic, B. et al. Fast image reconstruction with L2‐regularization. 
J Magn Reson Imaging 40, 181–191 (2014).](https://doi.org/10.1002/jmri.24365)  

**LSQR**  
[Li, W., Wu, B. & Liu, C. Quantitative susceptibility mapping of human brain reflects spatial 
variation in tissue composition. 
Neuroimage 55, 1645–1656 (2011).](https://doi.org/10.1016/j.neuroimage.2010.11.088)  

**Star-QSM**  
[Wei, H. et al. Streaking artifact reduction for quantitative susceptibility mapping of sources with 
large dynamic range. NMR Biomed 28, 1294–1303 (2015).](https://doi.org/10.1002/nbm.3383)  

[Wei, H. et al. Imaging whole-brain cytoarchitecture of mouse with MRI-based quantitative 
susceptibility mapping. 
Neuroimage 137, 107–115 (2016).](https://doi.org/10.1016/j.neuroimage.2016.05.033)  

[Wei, H. et al. Investigating magnetic susceptibility of human knee joint at 7 Tesla. 
Magn Reson Med 78, 1933–1943 (2017).](https://doi.org/10.1002/mrm.26596)  

**FANSI**  
[Milovic, C., Bilgic, B., Zhao, B., Acosta-Cabronero, J. & Tejos, C. Fast nonlinear susceptibility 
inversion with variational regularization. 
Magn Reson Med 80, 814–821 (2018).](https://doi.org/10.1002/mrm.27073)  

[Bilgic, B. et al. Fast quantitative susceptibility mapping with L1‐regularization and automatic 
parameter selection. Magn Reson Med 72, 1444–1459 (2014).](https://doi.org/10.1002/mrm.25029)  

[Bilgic, B., Chatnuntawech, I., Langkammer, C. & Setsompop, K. Sparse methods for Quantitative 
Susceptibility Mapping. in (eds. Papadakis, M., Goyal, V. K. & Van De Ville, D.) 9597, 959711 
(SPIE, 2015).](https://doi.org/10.1117/12.2188535)

**MEDI**  
[Liu, T. et al. Morphology enabled dipole inversion (MEDI) from a single-angle acquisition: 
Comparison with COSMOS in human brain imaging. 
Magn Reson Med 66, 777–783 (2011).](https://doi.org/10.1002/mrm.22816)  

[Liu, J. et al. Morphology enabled dipole inversion for quantitative susceptibility mapping using 
structural consistency between the magnitude image and the susceptibility map. 
Neuroimage 59, 2560–2568 (2012).](https://doi.org/10.1016/j.neuroimage.2011.08.082)  

[Liu, Z., Spincemaille, P., Yao, Y., Zhang, Y. & Wang, Y. MEDI+0: Morphology enabled dipole 
inversion with automatic uniform cerebrospinal fluid zero reference for quantitative susceptibility 
mapping. Magn Reson Med 79, 2795–2803 (2018).](https://doi.org/10.1002/mrm.26946)


----------------------------------------------------------------------------------------------------

## Standalone description
--------------------------------------------------

### QSMHub (One-stop QSM processing)
==================================================

#### I/O panel
--------------------------------------------------
- Input:  
	(Option 1) 	Directory contains all and only mGRE DICOM data, including both magnitude and phase
			images  
	(Option 2) 	Directory contains both magnitude and phase 4D-NIfTI data([x,y,slice,time]) and qsmhub
			header mat file (e.g. qsmhub_header.mat, optional), the magnitude data filename must
			contain the string 'magn' (e.g. 'magn.nii.gz') and the phase data filename must contain
			the string 'phase' (e.g. 'phase.nii.gz')  

- Output:
	Directory to store all qsm_hub output (default: '/input/directory/output')  
  
- FSL brain extraction (optional):
	Simple FSL's BET script, which can roughly extract a brain yet not too accurate.
	User is recommended to perform accurate brain extraction outside the qsm_hub app (e.g. using
	FSL's bet which allows fine tuning for more accurate result)  
  
- Brain mask file:
	NIfTI file that brain voxel > 0.
	QSM relies on brain mask thorough all processing steps.
	If the input directory also contains NIfTI file that contains string 'mask' in filename (e.g.
	'mask.nii.gz'), it will be read automatically so you might skip this steps  

*Example  
A standard input directory contains the following files:  
- magn.nii.gz			(4D real image)  
- phase.nii.gz			(4D real image)  
- mask.nii.gz 			(3D mask image)  
- qsmhub_header.mat		(.mat file)*  

####	Total field recovery and phase unwrapping panel
--------------------------------------------------
- Method:  
	1. [**Laplacian**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
		very reliable method for phase unwrapping yet the output values may not be accurate  
		
	2. [**Laplacian STI suite**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
		Laplacian unwrapping implementation from STI Suite v3.0  
		
	3. [**3D best path**](https://doi.org/10.1364/AO.46.006623)   
		very robust region growing method yet only works in the DCCN cluster 
		(recommended if you use this toolbox in the DCCN cluster)  
		
	4. **Region growing**  
		MEDI toolbox implementation, might not work well with DICOM phase data 
		
	5. [**Graphcut**](https://doi.org/10.1109/TMI.2014.2361764)  
		graph-cut algorithm (not optimised with this toolbox)  
		
- Bipolar readout eddy current correction:  
	enable to correct the phase inconsistence between odd and even echoes and gradient field
	contribution by eddy current effect due to bipolar readout  
  
- Exclude unreliable voxels, Threshold:  
	enable to exclude low SNR voxels that can create strong artefact in susceptibility map (you may
	check with 'squirrel_fieldmapSD.nii.gz' to adjust the threshold)  
	
*Output  
- squirrel_phase_EC.nii.gz 		(eddy current corrected phase, optional)  
- squirrel_magn_EC.nii.gz		(eddy current corrected magnitude, optional)  
- squirrel_totalField.nii.gz 	(unwrapped total (background+local) field, in Hz)  
- squirrel_fieldMapSD.nii.gz 	(normalised field map standaed deviation)*  

####	Background field removal panel
--------------------------------------------------
- Method:
	1. [**LBV**](https://doi.org/10.1002/nbm.3064)  
		Laplacian boundary value approach to removal background field 
		
	2. [**PDF**](https://doi.org/10.1002/nbm.1670)  
		Projection onto dipole field
		
	3. [**RESHARP**](https://doi.org/10.1002/mrm.24765)  
		Regularisation enabled SHARP
		
	4. [**SHARP**](https://doi.org/10.1016/j.neuroimage.2010.10.070)  
		Sophisticated harmonic artefact reduction for phase data
		
	5. [**VSHARP STI suite**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
		STI suite v3.0 variable-kernel SHARP 
		
	6. [**VSHARP**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
	
	7. [**iHARPERELLA**](https://doi.org/10.1002/nbm.3056)  
		(not optimised for this toolbox)  

- Refine local field by 4th order 3D polynomial fit
	enable to remove residual B1(+ & -) contribution in the local field  

*Output  
- squirrel_localField.nii.gz 		(background-field-removed local (or tissue) field, in Hz)  
- squirrel_mask_final.nii.gz		(final brain mask, might be eroded depended on background field
									removal algorithms and 'exclude unrelaiable voxels' threshold value)*  

####	QSM panel
--------------------------------------------------
- Method:
	1. [**TKD**](https://doi.org/10.1002/mrm.22334)  
		Thresholded k-space division  
		
	2. [**Closed-form solution**](https://doi.org/10.1002/jmri.24365)  
		Closed-form solution with L2-norm regularisation  
		
	3. [**STI suite iLSQR**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
		STI suite v3.0 implementation of iterative LSQR approach  
		
	4. [**iLSQR**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
	
	5. [**FANSI**](https://doi.org/10.1002/mrm.27073)  
		Fast algorithm for nonlinear susceptibility inversion   
		
	6. [**Star**](https://doi.org/10.1002/nbm.3383)  
		STI suite v3.0 Star-QSM   
		
	7. [**MEDI**](https://doi.org/10.1002/mrm.26946)  
		Morphology enabled dipole inversion (MEDI+0) 
  
*Output  
- squirrel_QSM.nii.gz 			(quantitative susceptibility map, in ppm)*

--------------------------------------------------
### Phase unwrapping
==================================================

####	I/O panel
--------------------------------------------------
- Input:  
	(Option 1) 	Directory contains all and only mGRE DICOM data, including both magnitude and phase
			images  
	(Option 2) 	Directory contains both magnitude and phase 4D-NIfTI data([x,y,slice,time]) and qsmhub
			header mat file (e.g. qsmhub_header.mat, optional), the magnitude data filename must
			contain the string 'magn' (e.g. 'magn.nii.gz') and the phase data filename must contain
			the string 'phase' (e.g. 'phase.nii.gz')  
  
- Output:  
	Directory to store all qsm_hub output (default: /input/dir/output/)  
  
- FSL brain extraction (optional):  
	Simple FSL's BET script, which can roughly extract a brain yet not too accurate.
	User is recommended to perform accurate brain extraction outside the qsm_hub app (e.g. using
	FSL's bet which allows fine tuning for more accurate result)  
  
- Brain mask file:  
	NIfTI file that brain voxel > 0.
	QSM relies on brain mask thorough all processing steps.
	If the input directory also contains NIfTI file that contains string 'mask' in filename (e.g.
	'mask.nii.gz'), it will be read automatically so you might skip this steps  

####	Total field recovery and phase unwrapping panel
--------------------------------------------------
- Method:
	1. [**Laplacian**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
		very reliable method for phase unwrapping yet the output values are not accurate
		
	2. [**Laplacian STI suite**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
		Laplacian unwrapping implementation from STI Suite v3.0  
		
	3. [**3D best path**](https://doi.org/10.1364/AO.46.006623)   
		very robust region growing method yet only works in the DCCN cluster 
		(recommended if you use this toolbox in the DCCN cluster)
		
	4. **Region growing**  
		MEDI toolbox implementation, might not work well with DICOM phase data  
  
	5. [**Graphcut**](https://doi.org/10.1109/TMI.2014.2361764)  
		graph-cut algorithm 
		(not optimised with this toolbox) 
		
- Bipolar readout eddy current correction:  
	enable to correct the phase inconsistence between odd and even echoes and gradient field
	contribution by eddy current effect due to bipolar readout  
	
- Exclude unreliable voxels, Threshold:  
	enable to exclude low SNR voxels that can create strong artefact in susceptibility map (you may
	check with 'squirrel_fieldmapSD.nii.gz' to adjust the threshold)    
	
*Output    
- squirrel_phase_EC.nii.gz 		(eddy current corrected phase, optional)  
- squirrel_magn_EC.nii.gz		(eddy current corrected magnitude, optional)  
- squirrel_totalField.nii.gz 	(unwrapped total (background+local) field, in Hz)  
- squirrel_fieldMapSD.nii.gz 	(normalised field map standaed deviation)*    

--------------------------------------------------
### Background field removal
==================================================

####	I/O panel
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
	
*Example  
A standard input directory contains the following files:  
- totalField.nii.gz		(unwrapped total field map)  
- fieldMapSD.nii.gz		(normalised field map standard deviation, optional)  
- mask.nii.gz 			(3D mask image)  
- qsmhub_header.mat 	(.mat file)*  

####	Background field removal panel
--------------------------------------------------
- Method:
	1. [**LBV**](https://doi.org/10.1002/nbm.3064)  
		Laplacian boundary value approach to removal background field  

	2. [**PDF**](https://doi.org/10.1002/nbm.1670)  
		Projection onto dipole field  

	3. [**RESHARP**](https://doi.org/10.1002/mrm.24765)  
		Regularisation enabled SHARP  

	4. [**SHARP**](https://doi.org/10.1016/j.neuroimage.2010.10.070)  
		Sophisticated harmonic artefact reduction for phase data  

	5. [**VSHARP STI suite**](https://doi.org/10.1016/j.neuroimage.2010.11.088)   
		STI suite v3.0 variable-kernel SHARP 

	6. [**VSHARP**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  

	7. [**iHARPERELLA**](https://doi.org/10.1002/nbm.3056)   
		(not optimised with this toolbox)  
		
- Refine local field by 4th order 3D polynomial fit
	enable to remove residual B1(+ & -) contribution in the local field
	
*Output  
- squirrel_localField.nii.gz 	(background-field-removed local (or tissue) field, in Hz)  
- squirrel_mask_final.nii.gz	(final brain mask, might be eroded depended on background field removal
								algorithms and 'exclude unrelaiable voxels' threshold value)*  

--------------------------------------------------
### QSM
==================================================

####	I/O planel
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
	
*Example  
A standard input directory contains the following files:  
- localField.nii.gz 	(3D local field map, in Hz)  
- magn.nii.gz			(4D real image used to cerate weighting map for some QSM algorithm, [x,y,slice,time],
						optional)  
- mask_final.nii.gz 	(3D mask image)  
- qsmhub_header.mat 	(.mat file)*  

####	QSM panel
--------------------------------------------------
- Method:
	1. [**TKD**](https://doi.org/10.1002/mrm.22334)  
		Thresholded k-space division

	2. [**Closed-form solution**](https://doi.org/10.1002/jmri.24365)  
		closed-form solution with L2-norm regularisation

	3. [**STI suite iLSQR**](https://doi.org/10.1016/j.neuroimage.2010.11.088)  
		STI suite v3.0 implementation of iterative LSQR approach

	4. [**iLSQR**](https://doi.org/10.1016/j.neuroimage.2010.11.088)

	5. [**FANSI**](https://doi.org/10.1002/mrm.27073)  
		Fast algorithm for nonlinear susceptibility inversion

	6. [**Star**](https://doi.org/10.1002/nbm.3383)  
		STI suite v3.0 Star-QSM (recommended)

	7. [**MEDI**](https://doi.org/10.1002/mrm.26946)  
		Morphology enabled dipole inversion (MEDI+0)   

*Output     
- squirrel_QSM.nii.gz 			(quantitative susceptibility map, in ppm)*  
