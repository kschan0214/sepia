fullName = mfilename('fullpath');

currDir = fileparts(fullName);

bkgRemovalFOLDER = [currDir filesep 'background_removal/'];
addpath(bkgRemovalFOLDER);
addpath([bkgRemovalFOLDER 'iHARPERELLA']);
addpath([bkgRemovalFOLDER 'LBV']);
addpath([bkgRemovalFOLDER 'PDF']);
addpath([bkgRemovalFOLDER 'SHARP']);
addpath([bkgRemovalFOLDER 'VSHARP']);
addpath([bkgRemovalFOLDER 'RESHARP']);
addpath([bkgRemovalFOLDER 'utils']);

phaseUnwrapFOLDER = [currDir filesep 'phase_unwrap/'];
addpath(phaseUnwrapFOLDER);
addpath([phaseUnwrapFOLDER 'unwrapGraphcut']);
addpath([phaseUnwrapFOLDER 'unwrapJena']);
addpath([phaseUnwrapFOLDER 'unwrapLaplacian']);
addpath([phaseUnwrapFOLDER 'unwrapRegionGrowing']);

qsmAlgorithmFOLDER = [currDir filesep 'qsm_algorithm/'];
addpath(qsmAlgorithmFOLDER);
addpath([qsmAlgorithmFOLDER 'closedFormL2']);
addpath([qsmAlgorithmFOLDER 'iLSQR']);
addpath([qsmAlgorithmFOLDER 'iLSQR' filesep 'STISuite']);
addpath([qsmAlgorithmFOLDER 'TKD']);
addpath([qsmAlgorithmFOLDER 'FANSI']);
addpath([qsmAlgorithmFOLDER 'SingleStep']);
addpath([qsmAlgorithmFOLDER 'utils']);

utilsFOLDER = [currDir filesep 'utils/'];
addpath(utilsFOLDER);
addpath([utilsFOLDER 'nifti/NIfTI_20140122/']);
addpath([utilsFOLDER 'nifti/utils/']);
addpath([utilsFOLDER 'ReadDICOMMEDI/']);