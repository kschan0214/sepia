bkgRemovalFOLDER = 'background_removal/';
addpath(bkgRemovalFOLDER);
addpath([bkgRemovalFOLDER 'iHARPERELLA']);
addpath([bkgRemovalFOLDER 'LBV']);
addpath([bkgRemovalFOLDER 'PDF']);
addpath([bkgRemovalFOLDER 'SHARP']);
addpath([bkgRemovalFOLDER 'VSHARP']);
addpath([bkgRemovalFOLDER 'RESHARP']);
addpath([bkgRemovalFOLDER 'utils']);

phaseUnwrapFOLDER = 'phase_unwrap/';
addpath(phaseUnwrapFOLDER);
addpath([phaseUnwrapFOLDER 'unwrapGraphcut']);
addpath([phaseUnwrapFOLDER 'unwrapJena']);
addpath([phaseUnwrapFOLDER 'unwrapLaplacian']);
addpath([phaseUnwrapFOLDER 'unwrapRegionGrowing']);

qsmAlgorithmFOLDER = 'qsm_algorithm/';
addpath(qsmAlgorithmFOLDER);
addpath([qsmAlgorithmFOLDER 'closedFormL2']);
addpath([qsmAlgorithmFOLDER 'iLSQR']);
addpath([qsmAlgorithmFOLDER 'iLSQR' filesep 'STISuite']);
addpath([qsmAlgorithmFOLDER 'TKD']);
addpath([qsmAlgorithmFOLDER 'FANSI']);
addpath([qsmAlgorithmFOLDER 'SingleStep']);
addpath([qsmAlgorithmFOLDER 'utils']);

addpath('utils');
addpath('utils/nifti/NIfTI_20140122/');
addpath('utils/nifti/utils/');
addpath('utils/ReadDICOMMEDI/');