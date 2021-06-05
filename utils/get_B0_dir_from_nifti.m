%% B0_dir = get_B0_dir_from_nifti(nii)
%
% Input
% --------------
% nii           : NIfTI format structure from load_untouch_nii.m
% B0_dir        : default B0 direction (optional, default in z-direction [0;0;1])
%
% Output
% --------------
% B0_dir        : main magnetic field direction with respect to imaging
%                 field of view
%
% Description: compute the B0 direction w.r.t. FOV
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 21 Jan 2021
% Date modified:
%
%
function B0_dir = get_B0_dir_from_nifti(nii,B0_dir)

if nargin < 2
    % assuming B0 is in the z-direction
    B0_dir = [0;0;1];
end

% first element of quaternion
a = sqrt(1 - nii.hdr.hist.quatern_b^2 - nii.hdr.hist.quatern_c^2 - nii.hdr.hist.quatern_d^2);

% transform from qiaternion representation to rotation matrix
rotmat_nifti = qGetR([a, nii.hdr.hist.quatern_b,nii.hdr.hist.quatern_c,nii.hdr.hist.quatern_d]);

% apply rotation matrix to B0
B0_dir = rotmat_nifti \ B0_dir;

end