function add_path_STIsuitev3(STISuite_dir)

addpath([STISuite_dir filesep 'Core_Functions_P']);
addpath([STISuite_dir filesep 'GUI_Functions_P']);
addpath([STISuite_dir filesep 'Support_Functions']);
addpath(genpath([STISuite_dir filesep 'Support_Functions' filesep 'qsm_kiwi_1']));
addpath(genpath([STISuite_dir filesep 'Support_Functions' filesep 'SpaRSA']));
addpath(genpath([STISuite_dir filesep 'Support_Functions' filesep 'wavelet_src']));

end