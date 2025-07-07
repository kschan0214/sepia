%% get_set_qsm_heidi(h,mode,input)
%
% Input
% --------------
% h             : structure contains all handles of SEPIA
% mode          : 'set' - extract information from GUI to config file
%                 'get' - extract information from config file to GUI
% input         : if mode is 'set' then input should be fid
%                 if mode is 'get' then input should be confige file text
%
% Description: Information communication between config file and GUI
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 June 2025
% Date last modified:
%
%
function get_set_qsm_heidi(h,mode,input)

str_pattern = {'.qsm.tolerance',...
               '.qsm.maxiter',...
               '.qsm.residualWeighting',...
               '.qsm.PostProcCone_threshold',...
               '.qsm.PostProcCone_tol',...
               '.qsm.PostProcCone_tolEnergy',...
               '.qsm.offsetUseBool',...
               '.qsm.isFourierDomainFormula',...
               '.qsm.TikhonovRegularizationSusceptibility',...
               '.qsm.solvingType',...
               '.qsm.DipoleFilter'};

action_handle = {h.qsm.HEIDI.edit.tol,...
                 h.qsm.HEIDI.edit.maxIter,...
                 h.qsm.HEIDI.edit.ResidualWeighting,...
                 h.qsm.HEIDI.edit.ProcProcConeThreshold,...
                 h.qsm.HEIDI.edit.ProcProcConeTol,...
                 h.qsm.HEIDI.edit.ProcProcConeTolEnergy,...
                 h.qsm.HEIDI.checkbox.isOffsetUse,...
                 h.qsm.HEIDI.checkbox.isFourierDomainFormula,...
                 h.qsm.HEIDI.popup.Tikhonov,...
                 h.qsm.HEIDI.popup.solver,...
                 h.qsm.HEIDI.popup.DipoleFilter};
           

switch lower(mode)
    case 'set'
        fid = input;
        
        for k = 1:6
            fprintf(fid,'algorParam%s = %s ;\n'         ,str_pattern{k},get(action_handle{k},	'String'));
        end

        for k = 7:8
            fprintf(fid,'algorParam%s = %i ;\n'             ,str_pattern{k},get(action_handle{k},	'Value'));
        end
        
        for k = 9:11
            fprintf(fid,'algorParam%s = ''%s'' ;\n'     ,str_pattern{k},action_handle{k}.String{action_handle{k}.Value,1});
        end
        
    case 'get'
        
        config_txt = input;
        
        % first 5 edit fields
        for k = 1:6
            pattern_curr    = str_pattern{k};
            val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
            set_non_nan_value(action_handle{k},'String',val)
        end

        % checkbox
        for k = 7:8
            pattern_curr    = str_pattern{k};
            val             = get_num_as_string(config_txt, pattern_curr, '=', ';');
            set_non_nan_value(action_handle{k}, 'Value', str2double(val))
        end

        % offsetUseBool
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch lower(val)
            case 'default'
                set_non_nan_value(action_handle{k},'Value',1)
            case 'partial gradient weighting'
                set_non_nan_value(action_handle{k},'Value',2)
            case 'laplacian'
                set_non_nan_value(action_handle{k},'Value',3)
        end

        % solver
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch lower(val)
            case 'default'
                set_non_nan_value(action_handle{k},'Value',1)
            case 'inversefiltering'
                set_non_nan_value(action_handle{k},'Value',2)
            case 'spatialdomaintv'
                set_non_nan_value(action_handle{k},'Value',3)
        end

        % dipole filter
        k = k+1;
        pattern_curr    = str_pattern{k};
        val             = get_string_as_string(config_txt, pattern_curr);
        switch lower(val)
            case 'default'
                set_non_nan_value(action_handle{k},'Value',1)
            case 'truncsingularvalues'
                set_non_nan_value(action_handle{k},'Value',2)
        end
        
end