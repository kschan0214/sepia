clear all;clc;close all;
% STEP 1: Import data
[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=...
Read_GE_DICOM('AXL_QSM');

%%%%% provide a Mask here if possible %%%%%%
if (~exist('Mask','var'))                     
    Mask = genMask(iField, voxel_size);
end

%%%%% provide a noise_level here if possible %%%%%%
if (~exist('noise_level','var'))
    noise_level = calfieldnoise(iField, Mask);
end

%%%%% normalize signal intensity by noise to get SNR %%%
iField = iField/noise_level;

%%%% Generate the Magnitude image %%%%
iMag = sqrt(sum(abs(iField).^2,4));

% STEP 2a: Field Map Estimation
%%%%%Estimate the frequency offset in each of the voxel using a 
%%%%%complex fitting %%%%
% if abs((TE(2)-TE(1))-(TE(3)-TE(2)))< 0.0002
%     [iFreq_raw N_std] = Fit_ppm_complex(iField);
% else
%     [iFreq_raw N_std] = Fit_ppm_complex_TE(iField);
% end
   
% STEP 2b: Spatial phase unwrapping %%%%
% iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
[water fat iFreq unwph_uf unwph N_std] = spurs_gc(iField,TE,CF,voxel_size);
%*** For high-resolution data with large matrix size, may try subsampling ***% 
% [water fat iFreq unwph_uf unwph N_std] = spurs_gc(iField,TE,CF,voxel_size,2);


% STEP 2c: Background Field Removal
%%%% Background field removal %%%%
[RDF shim] = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);


% STEP 3: Dipole Inversion
save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
     voxel_size delta_TE CF B0_dir;
%%%% run MEDI %%%%%
QSM = MEDI_L1('lambda',1000);

