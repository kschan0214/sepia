%% sepia_load_addons
%
% Description: a script provides add-on ability in SEPIA
% names
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 27 June 2020
% Date modified: 
%
% DO NOT change the variable name
% DO NOT change the order of the entities, add a new one at the end instead
%
%% find addons in these directories
addons_dir          = fullfile(SEPIA_HOME,'addons');
addons_unwrap_dir   = fullfile(addons_dir,'phase_unwrap');
addons_bfr_dir      = fullfile(addons_dir,'bfr');
addons_qsm_dir      = fullfile(addons_dir,'qsm');

%% Phase unwrapping addons
listing = dir(addons_unwrap_dir);

for klist = 3:length(listing)
    if listing(klist).isdir 
        curr_dir = fullfile(addons_unwrap_dir,listing(klist).name);
        if exist(fullfile(curr_dir,'addon_config.m'),'file')
            run(fullfile(curr_dir,'addon_config.m'))
            methodUnwrapName{end+1}        = addons.method;
            wrapper_Unwrap_function{end+1} = addons.wrapper_function;
            gui_unwrap_exclusion{end+1}    = addons.gui_exclude_voxel;
        end
    end

end

%% BFR addons
listing = dir(addons_bfr_dir);

for klist = 3:length(listing)
    if listing(klist).isdir 
        curr_dir = fullfile(addons_bfr_dir,listing(klist).name);
        if exist(fullfile(curr_dir,'addon_config.m'),'file')
            run(fullfile(curr_dir,'addon_config.m'))
            methodBFRName{end+1}        = addons.method;
            wrapper_BFR_function{end+1} = addons.wrapper_function;
            if ~isempty(addons.gui_method_panel)
                function_BFR_method_panel{end+1} = addons.gui_method_panel;
            end
            if ~isempty(addons.config_function)
                config_BFR_function{end+1} = addons.config_function;
            end
        end
    end

end

%% QSM addons
listing = dir(addons_qsm_dir);

for klist = 3:length(listing)
    if listing(klist).isdir 
        curr_dir = fullfile(addons_qsm_dir,listing(klist).name);
        if exist(fullfile(curr_dir,'addon_config.m'),'file')
            run(fullfile(curr_dir,'addon_config.m'))
            methodQSMName{end+1}        = addons.method;
            wrapper_QSM_function{end+1} = addons.wrapper_function;
            if ~isempty(addons.gui_method_panel)
                function_QSM_method_panel{end+1} = addons.gui_method_panel;
            end
            if ~isempty(addons.config_function)
                config_QSM_function{end+1} = addons.config_function;
            end
        end
    end

end
