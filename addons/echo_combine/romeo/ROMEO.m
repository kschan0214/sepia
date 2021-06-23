function totalField = ROMEO(phase, mag, mask, TE)
    tmp_dir = fullfile(tempdir, 'romeo_tmp'); % should create a suitable temporary directory on every machine
    mkdir(tmp_dir);
    % Input
    fn_phase = fullfile(tmp_dir, 'Phase.nii');
    fn_mag = fullfile(tmp_dir, 'Mag.nii'); 
    fn_mask = fullfile(tmp_dir, 'Mask.nii');
    % Output
    fn_unwrapped = fullfile(tmp_dir, 'Unwrapped.nii');
    fn_totalField = fullfile(tmp_dir, 'B0.nii');
    fn_phase_offset = fullfile(tmp_dir, 'phase_offset.nii');
    fn_corrected_phase = fullfile(tmp_dir, 'corrected_phase.nii');
    fn_romeo_settings = fullfile(tmp_dir, 'settings_romeo.txt');
    
    save_nii(make_nii(phase), fn_phase);
    save_nii(make_nii(mag), fn_mag);
    save_nii(make_nii(mask), fn_mask);
    path_to_binary = 'C:\Users\korbi\Desktop\romeo_win_3.2.0\bin';
    
    romeo_name = 'romeo';
    if ispc
        romeo_name = 'romeo.exe';
    end
    romeo_binary = fullfile(path_to_binary, romeo_name); 
    
    romeo_cmd = sprintf('%s %s -m %s -o %s -k %s -t %s -B --phase-offset-correction bipolar', romeo_binary, fn_phase, fn_mag, fn_unwrapped, fn_mask, mat2str(TE));
    success = system(romeo_cmd); % system command should work on every machine
    
    if success ~= 0
        error(['ROMEO unwrapping failed! Check input files for corruption in ' tmp_dir]);
    end
    
    totalField = load_nii_img_only(fn_totalField);
    delete(fn_totalField)
    delete(fn_phase)
    delete(fn_mag)
    delete(fn_mask)
    delete(fn_unwrapped)
    delete(fn_romeo_settings)
    delete(fn_phase_offset)
    delete(fn_corrected_phase)
    rmdir(tmp_dir)
end

