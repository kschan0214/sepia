%%%% Copyright WANGLAB 2011 All Rights Reserved. 
% created by Shuai 20110308 to read GE Dicom Files.
% for general GE DICOM images consisting of magnitude and phase
%
% Recommend: first go to the directory of case
% first : go to the directory of the case
% 
% DicomFloder: This folder should NOT contain anything else except the DICOM images belonging to a particular series that you want to process.
% [rawiField,voxel_size,matrix_size,CF,delta_TE,EchoN,Affine,B0_dir,minZ] =Read_GE_DICOM('SE3');
% based on Case data (complex data)
% Modified by Tian 20110519
% Modified by Shuai 20111020
function [real_p,voxel_size,matrix_size,CF,delta_TE,TE,NumEcho,Affine3D,B0_dir,minLoc]=Read_Bruker_DICOM_real(DicomFolder)
    
    filelist = dir(DicomFolder);
    i=1;
    while i<=length(filelist)
        if filelist(i).isdir==1
            filelist = filelist([1:i-1 i+1:end]);   % eliminate folders
        else
            i=i+1;
        end
    end

    fnTemp=[DicomFolder '/' filelist(1).name];
    
 
    
    info = dicominfo(fnTemp);
    matrix_size(1) = single(info.Width);
    matrix_size(2) = single(info.Height);
    matrix_size(3) = round(length(filelist)/info.ImagesInAcquisition);
    NumEcho = info.ImagesInAcquisition;
    %matrix_size(3) = single(info.Private_0021_104f);    
    %NumEcho = single(info.Private_0019_107e);
    voxel_size(1,1) = single(info.PixelSpacing(2));
    voxel_size(2,1) = single(info.PixelSpacing(1));
    voxel_size(3,1) = single(info.SliceThickness);
    CF = info.ImagingFrequency *1e6;
    TE = single(zeros([NumEcho 1]));
    real_p = single(zeros([matrix_size NumEcho]));
    minSlice = 1e10;
    maxSlice = -1e10;
    
    %for i = 1:length(filelist)/info.ImagesInAcquisition+1
    for i = 1:length(filelist)/info.ImagesInAcquisition
        info = dicominfo([DicomFolder '/' filelist(i).name]);
       
        if info.SliceLocation<minSlice
            minSlice = info.SliceLocation;
            minLoc = info.ImagePositionPatient;
        end
        if info.SliceLocation>maxSlice
            maxSlice = info.SliceLocation;
            maxLoc = info.ImagePositionPatient;
        end
    end
    
    Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
    Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
    B0_dir = Affine3D\[0 0 1]';
    

    for i = 1:length(filelist)
       
        info = dicominfo([DicomFolder '/' filelist(i).name]);
        
        slice = int32(round(sqrt(sum((info.ImagePositionPatient-minLoc).^2))/voxel_size(3)) +1);
        if TE(info.AcquisitionNumber)==0
            TE(info.AcquisitionNumber)=info.EchoTime*1e-3;
        end
        real_p(:,:,slice,info.AcquisitionNumber)  = single(dicomread([DicomFolder '/' filelist(i).name]))';%change 1&2 dimension ,magnitude
    end
    
    if length(TE)==1
        delta_TE = TE;
    else
        delta_TE = TE(2) - TE(1);
    end
    
end
