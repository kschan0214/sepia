function totalField = ROMEO_totalField(phase, mag, mask, parameters)
    % TODO set the path in SEPIA and retrieve it
    path_to_binary = 'C:\Users\korbi\Desktop\romeo_win_3.2.0\bin';

    % Should create a suitable temporary directory on every machine
    tmp_dir = fullfile(tempdir, 'romeo_tmp'); 
    mkdir(tmp_dir);
    
    % Input Files
    fn_phase = fullfile(tmp_dir, 'Phase.nii');
    fn_mag = fullfile(tmp_dir, 'Mag.nii'); 
    fn_mask = fullfile(tmp_dir, 'Mask.nii');
    save_nii(make_nii(phase), fn_phase);
    save_nii(make_nii(mag), fn_mag);
    save_nii(make_nii(mask), fn_mask);
    
    % Output Files
    fn_unwrapped = fullfile(tmp_dir, 'Unwrapped.nii');
    fn_totalField = fullfile(tmp_dir, 'B0.nii');
    
    % Always required parameters
    cmd_phase = [' -p ' fn_phase];
    cmd_output = [' -o ' fn_unwrapped];
    cmd_calculate_B0 = ' -B';
    cmd_mag = [' -m ' fn_mag];
    cmd_echo_times = [' -t ' mat2str(parameters.TE)];
    
    % Optional parameters
    cmd_mask = '';
    if strcmp(parameters.mask, 'SEPIA')
        cmd_mask = [' -k ' fn_mask];
    elseif strcmp(parameters.mask, 'nomask')
        cmd_mask = [' -k ' 'nomask'];
    elseif strcmp(paramaters.mask, 'romeomask')
        cmd_mask = [' -k ' 'robustmask'];
    end
    cmd_phase_offset_correction = '';
    if strcmp(parameters.phase_offset_correction, 'bipolar')
        cmd_phase_offset_correction = [' --phase-offset-correction ' 'bipolar'];
    elseif parameters.phase_offset_correction
        cmd_phase_offset_correction = [' --phase-offset-correction ' 'on'];
    end
    
    % Romeo binary name
    romeo_name = 'romeo';
    if ispc
        romeo_name = 'romeo.exe';
    end
    romeo_binary = fullfile(path_to_binary, romeo_name); 
    
    % Create romeo CMD command
    romeo_cmd = [romeo_binary cmd_phase cmd_mag cmd_mask cmd_output cmd_calculate_B0 cmd_echo_times cmd_phase_offset_correction];
    
    % Run romeo
    success = system(romeo_cmd); % system() call should work on every machine
    
    if success ~= 0
        error(['ROMEO unwrapping failed! Check input files for corruption in ' tmp_dir]);
    end
    
    % Load the calculated total field
    totalField = load_nii_img_only(fn_totalField);
    
    % Remove all temp output files and the temp folder
    rmdir(tmp_dir, 's')
end

