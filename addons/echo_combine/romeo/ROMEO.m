function [totalField, fieldmapUnwrapAllEchoes] = ROMEO(phase, mag, mask, parameters)
    % TODO set the path in SEPIA and retrieve it
    path_to_binary = 'C:\Users\korbi\Desktop\romeo_win_3.2.0\bin';

    % Should create a suitable temporary directory on every machine
    tmp_dir = fullfile(tempdir, 'romeo_tmp'); 
    mkdir(tmp_dir);
    
    % Input Files
    fn_phase = fullfile(tmp_dir, 'Phase.nii');
    fn_mask = fullfile(tmp_dir, 'Mask.nii');
    fn_mag = fullfile(tmp_dir, 'Mag.nii');
    
    phase_nii = make_nii(phase);
    phase_nii.hdr.dime.pixdim(2:4) = parameters.voxelSize;
    save_nii(phase_nii, fn_phase); % TODO set voxel size
    if parameters.useMag
        save_nii(make_nii(mag), fn_mag);
    end
    if strcmp(parameters.mask, 'file')    
        save_nii(make_nii(mask), fn_mask);
    end
    
    
    % Output Files
    fn_unwrapped = fullfile(tmp_dir, 'Unwrapped.nii');
    fn_totalField = fullfile(tmp_dir, 'B0.nii');
    
    % Always required parameters
    cmd_phase = [' -p ' fn_phase];
    cmd_output = [' -o ' fn_unwrapped];
    
    % Optional parameters
    cmd_calculate_B0 = '';
    if parameters.calculateB0
        cmd_calculate_B0 = ' -B';
    end
    cmd_mag = '';
    if parameters.useMag
        cmd_mag = [' -m ' fn_mag];
    end
    cmd_echo_times = [' -t ' mat2str(parameters.TE)];
    if strcmp(parameters.mask, 'file')
        cmd_mask = [' -k ' fn_mask];
    else
        cmd_mask = [' -k ' parameters.mask];
    end
    cmd_phase_offset_correction = [' --phase-offset-correction ' parameters.phaseOffsetCorrection];
    
    % Romeo binary name
    romeo_name = 'romeo';
    if ispc
        romeo_name = 'romeo.exe';
    end
    romeo_binary = fullfile(path_to_binary, romeo_name); 
    
    % Create romeo CMD command
    romeo_cmd = [romeo_binary cmd_phase cmd_mag cmd_output cmd_echo_times cmd_mask cmd_calculate_B0 cmd_phase_offset_correction];
    display(['ROMEO command: ' romeo_cmd]);
    
    % Run romeo
    success = system(romeo_cmd); % system() call should work on every machine
    
    if success ~= 0
        error(['ROMEO unwrapping failed! Check input files for corruption in ' tmp_dir]);
    end
    
    % Load the calculated output
    totalField = load_nii_img_only(fn_totalField);
    if parameters.isSaveUnwrappedEcho
        fieldmapUnwrapAllEchoes = load_nii_img_only(fn_unwrapped);
    else
        fieldmapUnwrapAllEchoes = [];
    end
    
    % Remove all temp output files and the temp folder
    rmdir(tmp_dir, 's')
end

