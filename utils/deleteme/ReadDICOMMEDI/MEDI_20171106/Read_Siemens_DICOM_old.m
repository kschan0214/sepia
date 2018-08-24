% for SIEMENS DICOM images consisting of magnitude and phase parts
%   [iField,voxel_size,matrix_size,CF,delta_TE,B0_dir] =Read_Siemens_DICOM(DicomFolder);
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
%   If there are open-ended fringe-lines on phase image,
%   1) save the uncombined data from each channel, put them into one folder
%   2) use this script to read the folder. 
%   The script will perform channel combination similarly to
%   MA. Bernstein et al. MRM 1994;32:330-334
%
%   Apdated from Tian Liu
%   Created by Shuai Wang on 2013.07.19 
%   Last modified by Tian Liu on 2014.02.13


function [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_Siemens_DICOM_old(DicomFolder)
    
    filelist = dir(DicomFolder);
    i=1;
    while i<=length(filelist)
        if filelist(i).isdir==1
            filelist = filelist([1:i-1 i+1:end]);   % skip folders
        else
            i=i+1;
        end
    end

    fnTemp=[DicomFolder '/' filelist(1).name];
    
 
    
    info = dicominfo(fnTemp);
    matrix_size(1) = single(info.Width);
    matrix_size(2) = single(info.Height);
    voxel_size(1,1) = single(info.PixelSpacing(1));
    voxel_size(2,1) = single(info.PixelSpacing(2));
    voxel_size(3,1) = single(info.SliceThickness);
    CF = info.ImagingFrequency *1e6;
    
    minSlice = 1e10;
    maxSlice = -1e10;
    NumEcho = 0;
    NumChannel = 0;
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
        if info.EchoNumber>NumEcho
            NumEcho = info.EchoNumber;
        end
        found = 0;
        for channel = 1:NumChannel
            if strcmp(info.Private_0051_100f, channelname{channel})
                found = 1;
                break;
            end
        end
        if found==0
            NumChannel = NumChannel + 1;
            channelname{NumChannel} = info.Private_0051_100f;
        end
    end
    matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3)) + 1;    
    Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
    Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
    B0_dir = Affine3D\[0 0 1]';

    iMag = single(zeros([matrix_size NumEcho NumChannel]));
    iPhase = single(zeros([matrix_size NumEcho NumChannel]));
    TE = single(zeros([NumEcho 1]));

    for i = 1:length(filelist)
        info = dicominfo([DicomFolder '/' filelist(i).name]);
        slice = int32(round(norm(info.ImagePositionPatient-minLoc)/voxel_size(3)) +1);
        if TE(info.EchoNumber)==0
            TE(info.EchoNumber)=info.EchoTime*1e-3;
        end
        
        for c = 1:NumChannel
            if strcmp(info.Private_0051_100f, channelname{c})
                channel = c;
            end
        end
        
        if (info.ImageType(18)=='P')|(info.ImageType(18)=='p')
            iPhase(:,:,slice,info.EchoNumber, channel)  = (single(dicomread([DicomFolder '/' filelist(i).name])')*info.RescaleSlope+info.RescaleIntercept)/single(info.LargestImagePixelValue)*pi;%phase
        elseif (info.ImageType(18)=='M')|(info.ImageType(18)=='m')
            iMag(:,:,slice,info.EchoNumber, channel)  = single(dicomread([DicomFolder '/' filelist(i).name])');%magnitude
        end
    end

    
    iField = iMag.*exp(-1i*iPhase);
    clear iMag iPhase;
    if NumChannel>1
        % combine multiple coils together, assuming the coil is the fifth dimension
        iField = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 NumEcho 1])),5);
        iField = sqrt(abs(iField)).*exp(1i*angle(iField));
    end

    
    if length(TE)==1
        delta_TE = TE;
    else
        delta_TE = TE(2) - TE(1);
    end
    
end