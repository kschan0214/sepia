%% function unwrappedField = UnwrapPhase_3DBestPath(wrappedField,mask,matrixSize)
%
% Input
% ----------------
%   wrappedField        : wrapped field map
%   mask                : brain mask
%   matrixSize          : images matrix size
%
% Output
% ----------------
%   unwrappedField      : unwrapped field
%
% Description: phase unwrapping using Jena's implementation of "Fast 
% three-dimensional phase-unwrapping algorithm based on sorting by 
% reliability following a noncontinuous path" 
% by Hussein Abdul-Rahman, Munther A. Gdeisat, David R. Burton, and 
% Michael J. Lalor, published in the Applied Optics, 
% Vol. 46, No. 26, pp. 6623-6635, 2007. (Only works in the DCCN cluster)
%
% Required dependency: Tools for NIfTI and ANALYZE image by Jimmy Shen
% (https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
% In future, this may be replaced by Matlab's niftiread/niftiwrite functions 
%
% This code is modified from the T2starAndFieldCalc.m from Jose P. Marques
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 29 June 2017
% Date last modified: 24 May 2018
%
function unwrappedField = UnwrapPhase_3DBestPath(wrappedField,mask,matrixSize)

temp = make_nii(single(wrappedField));
%     temp.img=temp.img;
% puts noise outside the mask to give freedom to the software in regions outside the mask;
temp.img = temp.img.*mask+(1-mask).*(rand(matrixSize(1:3))*2*pi-pi);
temp.hdr.dime.datatype=16;
temp.hdr.dime.bitpix=16;
save_nii(temp,'temp.hdr');
fn = mfilename('fullpath');
[pathstr,~,~] = fileparts(fn);
unix([pathstr '/JenaUnwrapDONDERS.sh temp.hdr']);
%     unix('sh  /home/rebelo/Documents/MATLAB/phase_unwrapping/test temp.hdr')
% unix('sh  /home/mrphys/kwocha/Tools/phase_unwrap/unwrapJena/unwrap temp.hdr')
unwrappedFieldnii = load_untouch_nii('uwtemp.hdr');
unwrappedField = unwrappedFieldnii.img;
system('rm temp.hdr temp.img uwtemp.img uwtemp.hdr qmtemp.img qmtemp.hdr temp.mat');
end