function unwrapped = ROMEO(phase, mag, mask)
    tmp_dir = fullfile(tempdir, 'romeo_tmp'); % should create a suitable temporary directory on every machine
    mkdir(tmp_dir);
    fn_phase = fullfile(tmp_dir, 'Phase.nii');
    fn_mag = fullfile(tmp_dir, 'Mag.nii'); 
    fn_mask = fullfile(tmp_dir, 'Mask.nii');
    fn_unwrapped = fullfile(tmp_dir, 'Unwrapped.nii');
    fn_romeo_settings = fullfile(tmp_dir, 'settings_romeo.txt');
    
    save_nii(make_nii(phase), fn_phase);
    save_nii(make_nii(mag), fn_mag);
    save_nii(make_nii(mask), fn_mask);
    path_to_binary = 'C:\Users\korbi\Desktop\SEPIA\romeo_win_3.1.4\romeo_win_3.1.4\bin';
    
    romeo_name = 'romeo';
    if ispc
        romeo_name = 'romeo.exe';
    end
    romeo_binary = fullfile(path_to_binary, romeo_name); 
    
    % TODO add mask (recompile romeo)
    romeo_cmd = sprintf('%s %s -m %s -o %s', romeo_binary, fn_phase, fn_mag, fn_unwrapped);
    system(romeo_cmd); % system command should work on every machine
    
    unwrapped = load_nii_img_only(fn_unwrapped);
    delete(fn_phase)
    delete(fn_mag)
    delete(fn_mask)
    delete(fn_unwrapped)
    delete(fn_romeo_settings)
    rmdir(tmp_dir)
end

