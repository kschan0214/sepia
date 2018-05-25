%% GetVentricleMaskIOWrapper(input, outputDir, maskFullName, TE)
%
% Input
% --------------
% input         : full path of mGRE magnitude NIfTI file
% outputDir     : output directory to ventricle mask
% maskFullName  : full path of mask NIfTI file
% TE            : Echo times, in second
%
% Description: Wrapper function to segment lateral ventricle
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 25 may 2018
% Date last modified:
%
%
function GetVentricleMaskIOWrapper(input, outputDir, maskFullName, TE)
prefix = 'squirrel_';

% add path
qsm_hub_AddMethodPath('medi_l1');

disp('Loading data...')

inputNifti = load_untouch_nii(input);

[~,~,voxelSize,~,~,~,~]=SyntheticQSMHubHeader(inputNifti);

magn = abs(inputNifti.img);

% load mask
if ~isempty(maskFullName)
% first read mask if file is provided
    mask = load_nii_img_only(maskFullName) > 0;
else 
	mask = max(magn,[],4)/max(magn(:)) > 0.05;
end

% R2* mapping
disp('Segmenting ventricle...');

r2s = arlo(TE,magn);
maskCSF = extract_CSF(r2s,mask,voxelSize)>0;

outputNiftiTemplate = inputNifti;
% make sure the class of output datatype is double
outputNiftiTemplate.hdr.dime.datatype = 4;
% remove the time dimension info
outputNiftiTemplate.hdr.dime.dim(5) = 1;

nii_maskCSF = make_nii_quick(outputNiftiTemplate,maskCSF);

save_untouch_nii(nii_maskCSF,[outputDir filesep prefix 'mask_ventricleCSF.nii.gz']);

disp('Done');

end

% handy function to save result to nifti format
function nii = make_nii_quick(template,img)
    nii = template;
    nii.img = img;
    nii.hdr.dime.datatype = 64;
end