%% function export_nii(img,fn,v_size)
%
% Description: Directly save to nii.gz format
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
function export_nii(img,fn,v_size)
load_module_NIfTI;

try filename=fn; catch; filename='Result';end
try voxel_size = v_size; catch; voxel_size=[1 1 1]; end

nii =make_nii(img,voxel_size);

save_nii(nii,[filename '.nii.gz']);
end
    