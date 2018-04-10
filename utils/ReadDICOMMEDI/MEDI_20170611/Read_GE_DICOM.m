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

function [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_GE_DICOM(DicomFolder)

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
isz=[matrix_size(1) matrix_size(2) matrix_size(3)*NumEcho];
iReal = single(zeros(isz));
iImag = single(zeros(isz));

TE = single(zeros([NumEcho 1]));
r_ImagePositionPatient = zeros(3, matrix_size(3)*NumEcho);
i_ImagePositionPatient = zeros(3, matrix_size(3)*NumEcho);
r_EchoNumber = zeros(matrix_size(3)*NumEcho,1);
minSlice = 1e10;
maxSlice = -1e10;


rctr = 0; ictr=0;
progress='';
for i = 1:length(filelist)
    filename=[DicomFolder '/' filelist(i).name];
    info2 = get_dcm_tags(filename);
    
    if info2.SliceLocation<minSlice
        minSlice = info2.SliceLocation;
        minLoc = info2.ImagePositionPatient;
    end
    if info2.SliceLocation>maxSlice
        maxSlice = info2.SliceLocation;
        maxLoc = info2.ImagePositionPatient;
    end
    if info2.EchoNumber>NumEcho
        TE = [TE zeros([info2.EchoNumber - NumEcho 1])];
        NumEcho = info2.EchoNumber;
    end
    if TE(info2.EchoNumber)==0
        TE(info2.EchoNumber)=info2.EchoTime*1e-3;
    end
    if mod(info2.InstanceNumber,3)==2
        rctr = rctr + 1;
        for ii=1:length(progress); fprintf('\b'); end
        progress=sprintf('Reading file %d', rctr+ictr);
        fprintf(progress);
        r_ImagePositionPatient(:,rctr) = info2.ImagePositionPatient;
        r_EchoNumber(rctr) = info2.EchoNumber;
        iReal(:,:,rctr)  = single(dicomread(filename));%magnitude
    elseif mod(info2.InstanceNumber,3)==0
        ictr = ictr + 1;
        for ii=1:length(progress); fprintf('\b'); end
        progress=sprintf('Reading file %d', rctr+ictr);
        fprintf(progress);
        i_ImagePositionPatient(:,ictr) = info2.ImagePositionPatient;
        i_EchoNumber(ictr) = info2.EchoNumber;
        iImag(:,:,ictr)  = single(dicomread(filename));%magnitude
    end
end
fprintf('\n');

matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3))+1 ;

Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
B0_dir = Affine3D\[0 0 1]';

sz=size(r_ImagePositionPatient); sz(1)=1;
minLoc=repmat(minLoc, sz);

r_slice = int32(round(sqrt(sum((r_ImagePositionPatient-minLoc).^2,1))/voxel_size(3)) +1);
r_ind = sub2ind([matrix_size(3) NumEcho], r_slice(:), int32(r_EchoNumber(:)));
iReal(:,:,r_ind)=iReal;
i_slice = int32(round(sqrt(sum((i_ImagePositionPatient-minLoc).^2,1))/voxel_size(3)) +1);
i_ind = sub2ind([matrix_size(3) NumEcho], i_slice(:), int32(i_EchoNumber(:)));
iImag(:,:,i_ind)=iImag;

%
%     for i = 1:length(filelist)
%         info = dicominfo([DicomFolder '/' filelist(i).name]);
%         slice = int32(round(sqrt(sum((info.ImagePositionPatient-minLoc).^2))/voxel_size(3)) +1);
%
%
%     end


iField = reshape(complex(iReal, iImag), ...
    [matrix_size(1) matrix_size(2) matrix_size(3) NumEcho]);
iField = permute(iField,[2 1 3 4 5]); %This is because the first dimension is row in DICOM but COLUMN in MATLAB
iField(:,:,1:2:end,:) = -iField(:,:,1:2:end,:);
if length(TE)==1
    delta_TE = TE;
else
    delta_TE = TE(2) - TE(1);
end

end

function info = get_dcm_tags(filename)
attrs=dicomattrs(filename);
for t={'SliceLocation','ImagePositionPatient',...
        'EchoTime','EchoNumber','InstanceNumber'}
    t=char(t);
    [gr, el] = dicomlookup(t);
    for i=1:length(attrs);
        if (attrs(i).Group==gr)&&(attrs(i).Element==el)
            eval(['info.' t '= sscanf(char(attrs(i).Data), ''%f\\'');'])
            break;
        end
    end
end
end