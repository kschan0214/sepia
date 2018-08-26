%write bruker ifield


function [iField, CF, B0_dir, Affine3D, TE, delta_TE, matrix_size, voxel_size] = Read_Bruker_DICOM(real_folder, imag_folder)

[real_p,voxel_size,matrix_size,CF,delta_TE,TE,NumEcho,Affine3D,B0_dir,minLoc]=Read_Bruker_DICOM_data(real_folder);
[imag_p,voxel_size,matrix_size,CF,delta_TE,TE,NumEcho,Affine3D,B0_dir,minLoc]=Read_Bruker_DICOM_data(imag_folder);
iField = real_p + 1i*imag_p;


% Write_iField(iField, 'iField');

B0_dir=round(B0_dir);
B0_dir=B0_dir./(B0_dir+eps);

% save iField.mat iField CF B0_dir Affine3D TE delta_TE matrix_size voxel_size

% save_path=pwd;
% save_parameter_txt(save_path);
% save('iField', 'iField', '-v7.3');

end





