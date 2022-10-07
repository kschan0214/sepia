%% [chi] = Wrapper_QSM_QSMnet(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% localField    : local field map (tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% chi           : magnetic susceptibility map, in ppm
%
% Description: This is a wrapper function to access NDI for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 October 2022
% Date modified: 
%
%
function [chi] = Wrapper_QSM_LPCNN(localField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
sepia_universal_variables;

% get algorithm parameters
% algorParam  = check_and_set_algorithm_default(algorParam);
% method      = algorParam.qsm.method;

% get extra data such as magnitude/weights/B0 direction/TE/etc.
headerAndExtraData = check_and_set_SEPIA_header_data(headerAndExtraData);
b0dir = headerAndExtraData.sepia_header.B0_dir;
b0    = headerAndExtraData.sepia_header.B0;

addon_lpcnn_dir = fileparts(mfilename('fullpath'));

% add path
sepia_addpath;

setup_LPCNN_environment;
inference_template  = fullfile(addon_lpcnn_dir, 'inference_template.py');
to_numpy_template   = fullfile(addon_lpcnn_dir, 'to_numpy_template.py');


temp_dir        = fullfile(pwd,'temp_lpcnn_dir');
temp_output_dir = fullfile(pwd, 'LPCNN','test_result','lpcnn_test_Bmodel');
to_numpy_fn     = fullfile(temp_dir, 'to_numpy.py'); 
inference_fn  	= fullfile(temp_dir, 'inference.py'); 

dipolar_mat_fn  = fullfile(temp_dir, 'dipole.mat');
dipolar_npy_fn  = fullfile(temp_dir, 'dipole.npy');
phase_fn        = fullfile(temp_dir, 'phase.nii.gz');
mask_fn         = fullfile(temp_dir, 'mask.nii.gz');

phase_list      = fullfile(temp_dir, 'phase.txt');
mask_list       = fullfile(temp_dir, 'mask.txt');
dipole_list     = fullfile(temp_dir, 'dipole.txt');

mkdir(temp_dir)
mkdir(temp_output_dir)
%% preparation
% clear GPU
gpuDevice([]);

% LPCNN only works with matrix with same in plane size
dims = size(localField);
if dims(2) > dims(1)
    zeropad_line = (dims(2) - dims(1))/2;
    localField_pad = zeros(dims(2),dims(2),dims(3));
    localField_pad(1+zeropad_line:end-zeropad_line,:,:) = localField;
elseif dims(2) < dims(1)
    zeropad_line = (dims(1) - dims(2))/2;
    localField_pad = zeros(dims(1),dims(1),dims(3));
    localField_pad(:,1+zeropad_line:end-zeropad_line,:) = localField;
end

mask_pad = localField_pad ~=0;
[C,~] = DipoleKernel(size(localField_pad),voxelSize,b0dir);

% export files for LPCNN
save(dipolar_mat_fn,'C')
export_nii(localField_pad,phase_fn);
export_nii(double(mask_pad),mask_fn);

%% main
% convert dipole from .mat to .npy file
copyfile(to_numpy_template,to_numpy_fn);
% should work for at least unix system
system(sprintf('sed -i "16i input_fn = ''%s''" %s',     dipolar_mat_fn,to_numpy_fn));
system(sprintf('sed -i "17i output_fn = ''%s''" %s',    dipolar_npy_fn,to_numpy_fn));

% run python
to_numpy_cmd = [python_interpreter ' ' to_numpy_fn];
system(to_numpy_cmd);

% prepare file lists
fid = fopen(phase_list,'w');
fprintf(fid,'%s\n',phase_fn);
fclose(fid);

fid = fopen(mask_list,'w');
fprintf(fid,'%s\n',mask_fn);
fclose(fid);

fid = fopen(dipole_list,'w');
fprintf(fid,'%s\n',dipolar_npy_fn);
fclose(fid);

% prepare inference.py
copyfile(inference_template,inference_fn);
% should work for at least unix system
system(sprintf('sed -i "16i LPCNN_HOME = Path(''%s'')" %s',	LPCNN_HOME,inference_fn));

% run LPCNN
lpcnn_cmd = [python_interpreter ' ' inference_fn ' --tesla ' num2str(b0) ' --number 1 --phase_file ' phase_list ' --dipole_file ' dipole_list ' --mask_file ' mask_list ' --resume_file ' checkpoint_fn];
system(lpcnn_cmd);

% load output back to matlab
chi = load_nii_img_only(fullfile(temp_output_dir, 'lpcnn_example_qsm.nii.gz'));

% post crop
if dims(2) > dims(1)
    chi = chi(1+zeropad_line:end-zeropad_line,:,:);
    
elseif dims(2) < dims(1)
    chi = chi(:,1+zeropad_line:end-zeropad_line,:);
end

chi = chi .* mask;

% remove all intermediate folders
rmdir(temp_dir, 's')
rmdir(fullfile(pwd, 'LPCNN'), 's')

end

% %% set default parameter if not specified
% function algorParam2 = check_and_set_algorithm_default(algorParam)
% 
% algorParam2 = algorParam;
% % 
% % try algorParam2.qsm.tol         = algorParam.qsm.tol;       catch; algorParam2.qsm.tol      = 1; end
% % try algorParam2.qsm.stepSize    = algorParam.qsm.stepSize;  catch; algorParam2.qsm.stepSize = 1; end
% % try algorParam2.qsm.maxiter     = algorParam.qsm.maxiter;   catch; algorParam2.qsm.maxiter  = 200; end
% 
% end