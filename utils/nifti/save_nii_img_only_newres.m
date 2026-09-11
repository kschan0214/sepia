%% function save_nii_img_only_newres(headerfilename,savefilename,img,datatype,voxelSize)
%
% Description: Wrapper script for saving with target image header
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 January 2017
% Date modified: 28 June 2020 (v0.8.0)
%
function save_nii_img_only_newres(headerfilename,savefilename,img,datatype,voxelSize)
% savefilename_o = savefilename;
% if verLessThan('matlab','9.3')
    if nargin < 4 || isempty(datatype); datatype = 16;   end
    if nargin < 5; voxelSize = [];  end

    nii = load_untouch_nii(headerfilename);

    % update image size
    nii.hdr.dime.dim(2:4) = size(img,1:3);
    % update voxel size if needed
    if ~isempty(voxelSize); nii.hdr.dime.pixdim(2:4) = voxelSize; end

    nii.img = single(img);
    nii.hdr.dime.datatype = datatype;
    nii.hdr.dime.dim(5) = size(img,4);
    nii.hdr.dime.dim(1) = ndims(img);

    nii.hdr.dime.scl_inter = 0;
    nii.hdr.dime.scl_slope = 1;

    save_untouch_nii(nii,savefilename);
end