function [phase,voxel_size,matrix_size,CF,delta_TE,TE,NumEcho,Affine3D,B0_dir,minLoc]=Read_Bruker_DICOM_data(realDir)

[real_p,voxel_size,matrix_size,CF,delta_TE,TE,NumEcho,Affine3D,B0_dir,minLoc]=Read_Bruker_DICOM_real(realDir);
real_p = real_p*((pi)/max(real_p(:)));

% [imag_p,voxel_size,matrix_size,CF,delta_TE,TE,NumEcho,Affine3D,B0_dir,minLoc]=Read_Bruker_DICOM_imag(imagDir);

% chopping (as per Pascal):
% for E = 1:length(TE)
%     for s = 1:matrix_size(3)
%         p=squeeze(real_p(:,:,s,E));
% 
%         p(1:2:end)=-p(1:2:end);
%         p(:,1:2:end)=-p(:,1:2:end);
%         real_p(:,:,s,E) = p;
%         
%     end
% end

%fftshift instead
if length(TE) > 1
    real_p = ifft(fftshift(fft(real_p, [], 1), 1), [], 1);
    real_p = ifft(fftshift(fft(real_p, [], 2), 2), [], 2);
    real_p = ifft(fftshift(fft(real_p, [], 3), 3), [], 3);
end    
    phase = real_p;

% 
% 
% % imag_p=double(imag_p);
% % % imag_p=squeeze(imag_p(:,:,zind));
% % imag_p(1:2:end,:)=-imag_p(1:2:end,:);
% % imag_p(:,1:2:end)=-imag_p(:,1:2:end);
% imag_p = ifft(fftshift(fft(imag_p, [], 1), 1), [], 1);
% imag_p = ifft(fftshift(fft(imag_p, [], 2), 2), [], 2);
% imag_p = ifft(fftshift(fft(imag_p, [], 3), 3), [], 3);

% rawiField=real_p+1i*imag_p;

%     for echo = 1:NumEcho
%          rawiField(:,:,:,echo) = ifft( fftshift( fft(rawiField(:,:,:,echo),[],3),3),[],3);
%     end
