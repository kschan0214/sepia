function [unwrapped, B0] = ROMEO(phase, parameters)
    romeo_path = 'C:\Users\korbi\Desktop\romeo_win_3.2.0\bin';
    romeo_name = 'romeo';
    if ispc
        romeo_name = 'romeo.exe';
    end
    romeo_binary = fullfile(romeo_path, romeo_name); 
    
    output_dir = parameters.output_dir;
    
    % Input Files
    fn_phase = fullfile(output_dir, 'Phase.nii');
    fn_mask = fullfile(output_dir, 'Mask.nii');
    fn_mag = fullfile(output_dir, 'Mag.nii');
    
    phase_nii = make_nii(phase);
    if isfield(parameters, 'voxel_size')
        phase_nii.hdr.dime.pixdim(2:4) = parameters.voxel_size;
    end
    save_nii(phase_nii, fn_phase);
    if ~isempty(parameters.mag)
        save_nii(make_nii(parameters.mag), fn_mag);
    end
    if isnumeric(parameters.mask)
        save_nii(make_nii(parameters.mask), fn_mask);
    end
    
    
    % Output Files
    fn_unwrapped = fullfile(output_dir, 'Unwrapped.nii');
    fn_total_field = fullfile(output_dir, 'B0.nii');
    
    % Always required parameters
    cmd_phase = [' -p ' fn_phase];
    cmd_output = [' -o ' fn_unwrapped];
    
    % Optional parameters
    cmd_calculate_B0 = '';
    if parameters.calculate_B0
        cmd_calculate_B0 = ' -B';
    end
    cmd_mag = '';
    if ~isempty(parameters.mag)
        cmd_mag = [' -m ' fn_mag];
    end
    cmd_echo_times = [' -t ' mat2str(parameters.TE)];
    if isnumeric(parameters.mask)
        cmd_mask = [' -k ' fn_mask];
    else
        cmd_mask = [' -k ' parameters.mask];
    end
    cmd_phase_offset_correction = [' --phase-offset-correction ' parameters.phase_offset_correction];
    additional_flags = '';
    if isfield(parameters, 'additional_flags')
        additional_flags = parameters.additional_flags;
    end
    
    % Create romeo CMD command
    romeo_cmd = [romeo_binary cmd_phase cmd_mag cmd_output cmd_echo_times cmd_mask cmd_calculate_B0 cmd_phase_offset_correction additional_flags];
    disp(['ROMEO command: ' romeo_cmd])
    
    % Run romeo
    success = system(romeo_cmd); % system() call should work on every machine
    
    if success ~= 0
        error(['ROMEO unwrapping failed! Check input files for corruption in ' output_dir]);
    end
    
    % Load the calculated output
    B0 = [];
    unwrapped = [];
    if parameters.calculate_B0
        B0 = load_nii_img_only(fn_total_field);
    end
    if ~parameters.no_unwrapped_output
        unwrapped = load_nii_img_only(fn_unwrapped);
    end
end

