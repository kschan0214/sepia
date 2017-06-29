%% function unwrappedField = unwrapJena(wrappedField,mask,matrixSize)
%
% Description: phase unwrapping using Jena's library (Only works in the cluster)
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
% This code is modified from the T2starAndFieldCalc.m from Jose P. Marques
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 29 June 2017
% Date last modified:
%
function unwrappedField = unwrapJena(wrappedField,mask,matrixSize)
load_module_NIFTI;

temp = make_nii(single(wrappedField));
%     temp.img=temp.img;
% puts noise outside the mask to give freedom to the software in regions outside the mask;
temp.img = temp.img.*mask+(1-mask).*(rand(matrixSize(1:3))*2*pi-pi);
temp.hdr.dime.datatype=16;
temp.hdr.dime.bitpix=16;
save_nii(temp,'temp.hdr');
%     unix('sh  /home/rebelo/Documents/MATLAB/phase_unwrapping/test temp.hdr')
unix('sh  /home/mrphys/kwocha/Tools/phase_unwrap/unwrapJena/unwrap temp.hdr')
unwrappedField = load_nii('uwtemp.hdr');
end