%% function save_nii_img_only(headerfilename,savefilename,images)
%
% Description: Wrapper script for saving with target image header
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 January 2017
% Date modified: 28 June 2020 (v0.8.0)
%
function save_nii_img_only(headerfilename,savefilename,img,datatype)
% savefilename_o = savefilename;
% if verLessThan('matlab','9.3')
    if nargin < 4
        datatype = 16;
    end

    nii = load_untouch_nii(headerfilename);
    nii.img = single(img);
    nii.hdr.dime.datatype = datatype;
    nii.hdr.dime.dim(5) = size(img,4);
    nii.hdr.dime.dim(1) = ndims(img);
    % if size(img,4) > 1
    %     nii.hdr.dime.dim(1) = 4;
    % else
    %     nii.hdr.dime.dim(1) = 3;
    % end

    % assume the input image contains the true values
    nii.hdr.dime.scl_inter = 0;
    nii.hdr.dime.scl_slope = 1;

    % if nii.hdr.dime.datatype ==4
    %     nii.hdr.dime.datatype = 16;
    % end
    % nii.img = images;

    save_untouch_nii(nii,savefilename);
% else
%     try 
%         info = niftiinfo(headerfilename);
% 
%         if nargin < 4
%             info.Datatype = 'single';
%             img = single(img);
%         end
%         [filepath,filename,ext] = fileparts(savefilename);
%         savefilename = [filepath filesep filename];
%         info.ImageSize = size(img);
%         info.PixelDimensions = info.PixelDimensions(1:ndims(img));
%         info.raw.dim(1) = ndims(img);
%         info.raw.dim(5) = size(img,4);
%         pixdim = zeros(size(info.raw.pixdim));
%         pixdim(2:length(info.PixelDimensions)+1) = info.raw.pixdim(2:length(info.PixelDimensions)+1);
%         info.raw.pixdim = pixdim;
%         niftiwrite(img,savefilename,info,'Compressed',true);
%         
%         if nargin < 4
%         datatype = 16;
%         end
%     catch
%         savefilename = savefilename_o;
%         if nargin < 4
%             datatype = 16;
%         end
%         nii = load_untouch_nii(headerfilename);
%         nii.img = single(img);
%         nii.hdr.dime.datatype = datatype;
%         nii.hdr.dime.dim(5) = size(img,4);
%         nii.hdr.dime.dim(1) = ndims(img);
%         % if size(img,4) > 1
%         %     nii.hdr.dime.dim(1) = 4;
%         % else
%         %     nii.hdr.dime.dim(1) = 3;
%         % end
% 
%         % assume the input image contains the true values
%         nii.hdr.dime.scl_inter = 0;
%         nii.hdr.dime.scl_slope = 1;
% 
%         % if nii.hdr.dime.datatype ==4
%         %     nii.hdr.dime.datatype = 16;
%         % end
%         % nii.img = images;
% 
%         save_untouch_nii(nii,savefilename);
%     end
% end
end