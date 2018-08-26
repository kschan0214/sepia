function [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_DICOM(DicomFolder)
    
    create_dicomattrs;

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
    manufacturer = info.Manufacturer;
    institute = info.InstitutionName;
    
    
    
    switch manufacturer
    
        case 'SIEMENS'
            [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_Siemens_DICOM(DicomFolder);

            disp('SIEMENS READ');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PART%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'GE MEDICAL SYSTEMS'
    
            if ~isempty(strfind(lower(institute), 'tongji')) || ~isempty(strfind(lower(institute), 'zju')) || ~isempty(strfind(lower(institute), 'aff2'))
                
                matrix_size(1) = single(info.Width);
                matrix_size(2) = single(info.Height);
                tmp = single(info.Private_0021_104f);
                matrix_size(3) = tmp(1);    
                NumEcho = single(info.Private_0019_107e);
                voxel_size(1,1) = single(info.PixelSpacing(1));
                voxel_size(2,1) = single(info.PixelSpacing(2));
                voxel_size(3,1) = single(info.SpacingBetweenSlices);
                CF = info.ImagingFrequency *1e6;
                iMag = single(zeros([matrix_size NumEcho]));
                iPhase = single(zeros([matrix_size NumEcho]));
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

                Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
                Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
                B0_dir = Affine3D\[0 0 1]';
                
                for i = 1:length(filelist)
                    info = dicominfo([DicomFolder '/' filelist(i).name]);
                    slice = int32(round(sqrt(sum((info.ImagePositionPatient-minLoc).^2))/voxel_size(3)) +1);
                    if TE(info.EchoNumber)==0
                        TE(info.EchoNumber)=info.EchoTime*1e-3;
                    end
                    if mod(info.InstanceNumber,2)==1
                        iMag(:,:,slice,info.EchoNumber)  = single(dicomread([DicomFolder '/' filelist(i).name])');%magnitude
                    elseif mod(info.InstanceNumber,2)==0
                        largestImagePixelValue = single(info.LargestImagePixelValue);
                        iPhase(:,:,slice,info.EchoNumber)  = single(dicomread([DicomFolder '/' filelist(i).name])')/largestImagePixelValue*pi;%phase
                    end
                end


                iField = iMag .* exp(sqrt(-1)*iPhase);
                clear iMag iPhase;
                for echo = 1:NumEcho
                     iField(:,:,:,echo) = ifft( fftshift( fft(iField(:,:,:,echo),[],3),3),[],3);
                end
                if length(TE)==1
                    delta_TE = TE;
                else
                    delta_TE = TE(2) - TE(1);
                end


                disp('GE READ other');                
                
            else


                [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_GE_DICOM(DicomFolder);

                disp('GE READ');
                
            end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PHILIPS NEED
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MOD
        case 'Philips Medical Systems'
    
            matrix_size(1) = single(info.Width);
            matrix_size(2) = single(info.Height);
            voxel_size(1,1) = single(info.PixelSpacing(1));
            voxel_size(2,1) = single(info.PixelSpacing(2));
            voxel_size(3,1) = single(info.SliceThickness);
            CF = info.ImagingFrequency *1e6;

            minSlice = 1e10;
            maxSlice = -1e10;
            NumEcho = 0;
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
            end
            matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3)) + 1;    
            Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
            Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
            B0_dir = Affine3D\[0 0 1]';

            iMag = single(zeros([matrix_size NumEcho]));
            iPhase = single(zeros([matrix_size NumEcho]));
            TE = single(zeros([NumEcho 1]));

            for i = 1:length(filelist)
                info = dicominfo([DicomFolder '/' filelist(i).name]);
                slice = int32(round(norm(info.ImagePositionPatient-minLoc)/voxel_size(3)) +1);
                if TE(info.EchoNumber)==0
                    TE(info.EchoNumber)=info.EchoTime*1e-3;
                end
                if (info.ImageType(18)=='P')|(info.ImageType(18)=='p')
                    iPhase(:,:,slice,info.EchoNumber)  = 1e-3*(single(dicomread([DicomFolder '/' filelist(i).name])')*info.RealWorldValueMappingSequence.Item_1.RealWorldValueSlope+info.RealWorldValueMappingSequence.Item_1.RealWorldValueIntercept);%phase
                elseif (info.ImageType(18)=='M')|(info.ImageType(18)=='m')
                    iMag(:,:,slice,info.EchoNumber)  = single(dicomread([DicomFolder '/' filelist(i).name])');%magnitude
                end
            end


            iField = iMag.*exp(-1i*iPhase);
            clear iMag iPhase;
            if length(TE)==1
                delta_TE = TE;
            else
                delta_TE = TE(2) - TE(1);
            end

            disp('PHILIPS READ');
        otherwise
            
            disp('LOADING FAILED');





       end
            
            
            
end