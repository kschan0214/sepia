% for GE SWAN images consisting of magnitude, real and imaginary parts
%   [iField,voxel_size,matrix_size,CF,delta_TE,B0_dir] =Read_GE_DICOM(DicomFolder);
%
%   output
%   iField - the multi-echo complex MR image
%   voxel_size - size of the voxel in mm
%   matrix_size - dimension of the field of view
%   CF - central frequency in Hz
%   delta_TE - TE spacing in s
%   TE - echo times in s
%   B0_dir - direction of the B0 field
%
%   input
%   DicomFloder: This folder should NOT contain anything else 
%                except the DICOM images belonging to a particular
%                series that you want to process.
%
%   Apdated from Tian Liu
%   Created by Shuai Wang on 2011.03.08 
%   Modified by Tian Liu on 2011.05.19
%   Last modified by Tian Liu on 2013.07.24

function [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_GE_SWAN(DicomFolder)
    
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
    matrix_size(3) = single(info.Private_0021_104f);    
    NumEcho = single(info.Private_0019_107e);
    voxel_size(1,1) = single(info.PixelSpacing(1));
    voxel_size(2,1) = single(info.PixelSpacing(2));
    voxel_size(3,1) = single(info.SpacingBetweenSlices);
    CF = info.ImagingFrequency *1e6;
    iReal = single(zeros([matrix_size NumEcho]));
    iImag = single(zeros([matrix_size NumEcho]));
    TE = single(zeros([NumEcho 1]));
    
    minSlice = 1e10;
    maxSlice = -1e10;

    for i = 1:length(filelist)
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

matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3))+1 ;    

    Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
    Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
    B0_dir = Affine3D\[0 0 1]';

    for i = 1:length(filelist)
        info = dicominfo([DicomFolder '/' filelist(i).name]);
        slice = int32(round(sqrt(sum((info.ImagePositionPatient-minLoc).^2))/voxel_size(3)) +1);
        if TE(info.EchoNumber)==0
            TE(info.EchoNumber)=info.EchoTime*1e-3;
        end
        if mod(info.InstanceNumber,3)==2
            iReal(:,:,slice,info.EchoNumber)  = single(dicomread([DicomFolder '/' filelist(i).name])');%magnitude
        elseif mod(info.InstanceNumber,3)==0
            iImag(:,:,slice,info.EchoNumber)  = single(dicomread([DicomFolder '/' filelist(i).name])');%magnitude
        end
    end

    
    iField = iReal + 1i*iImag;
    clear iReal iImag;
%     for echo = 1:NumEcho
%          iField(:,:,:,echo) = ifft( fftshift( fft(iField(:,:,:,echo),[],3),3),[],3);
%     end
    iField(:,:,1:2:end,:) = -iField(:,:,1:2:end,:);
    if length(TE)==1
        delta_TE = TE;
    else
        delta_TE = TE(2) - TE(1);
    end
    
end