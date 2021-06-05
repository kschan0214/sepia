
function set_config_Callback(config_filename,h)

sepia_universal_variables;

config_txt = fileread(config_filename);

%% to determine what is going to be modified
isModifyGeneral = false;
isModifyUnwrap  = false;
isModifyBFR     = false;
isModifyQSM     = false;

str_pattern	= '.general.';
idx         = regexp(config_txt,str_pattern,'once');
if ~isempty(idx)
    isModifyGeneral = true;
end
str_pattern	= '.unwrap.';
idx         = regexp(config_txt,str_pattern,'once');
if ~isempty(idx)
    isModifyUnwrap = true;
end
str_pattern	= '.bfr.';
idx         = regexp(config_txt,str_pattern,'once');
if ~isempty(idx)
    isModifyBFR = true;
end
str_pattern	= '.qsm.';
idx         = regexp(config_txt,str_pattern,'once');
if ~isempty(idx)
    isModifyQSM = true;
end 

%% Panel: I/O
if isModifyGeneral
    
% invert phase data
str_pattern	= '.general.isInvert';
val      	= get_num_as_string(config_txt, str_pattern,'=',';');
set_non_nan_value(h.dataIO.checkbox.invertPhase, 'Value', str2double(val));

% BET
str_pattern	= '.general.isBET';
val      	= get_num_as_string(config_txt, str_pattern,'=',';');
set_non_nan_value(h.dataIO.checkbox.brainExtraction, 'Value', str2double(val));
% trigger checkout callback
feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
if str2double(val)
    % -f
    str_pattern	= '.general.fractional_threshold';
    val      	= get_num_as_string(config_txt, str_pattern,'=',';');
    set_non_nan_value(h.dataIO.edit.fractionalThres, 'String', val);
    
    % -g
    str_pattern	= '.general.gradient_threshold';
    val      	= get_num_as_string(config_txt, str_pattern,'=',';');
    set_non_nan_value(h.dataIO.edit.gradientThres, 'String', val);
end

end

%% Panel: Total field
if isModifyUnwrap
    
% echo combination
str_pattern	= '.unwrap.echoCombMethod';
val	= get_string_as_string(config_txt, str_pattern);
if ~isnan(val)
    % matching algorithm name
    for k = 1:length(methodEchoCombineName)
        if strcmpi(val,methodEchoCombineName{k})
            val = k;
            break
        end
    end
    
    % change method popup manu
    set_non_nan_value(h.phaseUnwrap.popup.phaseCombMethod,'Value',val);
    % trigger popup callback to switch method panel
    feval(h.phaseUnwrap.popup.phaseCombMethod.Callback{1},h.phaseUnwrap.popup.phaseCombMethod,[],h)
end

% phase unwrap
str_pattern	= '.unwrap.unwrapMethod';
val	= get_string_as_string(config_txt, str_pattern);
if ~isnan(val)
    % matching algorithm name
    for k = 1:length(methodUnwrapName)
        if strcmpi(val,methodUnwrapName{k})
            val = k;
            break
        end
    end
    
    % change method popup manu
    set_non_nan_value(h.phaseUnwrap.popup.phaseUnwrap,'Value',val);
    % trigger popup callback to switch method panel
    feval(h.phaseUnwrap.popup.phaseUnwrap.Callback{1},h.phaseUnwrap.popup.phaseUnwrap,[],h)
end

% eddy corr
str_pattern	= '.unwrap.isEddyCorrect';
val      	= get_num_as_string(config_txt, str_pattern,'=',';');
set_non_nan_value(h.phaseUnwrap.checkbox.eddyCorrect, 'Value', str2double(val));

% exlusion threshold
str_pattern	= '.unwrap.excludeMaskThreshold';
val      	= get_num_as_string(config_txt, str_pattern,'=',';');
if isnan(val)
    % trigger checkbox
    set(h.phaseUnwrap.checkbox.excludeMask,'Value',0);
    feval(h.phaseUnwrap.checkbox.excludeMask.Callback{1},h.phaseUnwrap.checkbox.excludeMask,[],{h.phaseUnwrap.edit.excludeMask,h.phaseUnwrap.popup.excludeMethod},1);
else
    % trigger checkbox
    set(h.phaseUnwrap.checkbox.excludeMask,'Value',1);
    feval(h.phaseUnwrap.checkbox.excludeMask.Callback{1},h.phaseUnwrap.checkbox.excludeMask,[],{h.phaseUnwrap.edit.excludeMask,h.phaseUnwrap.popup.excludeMethod},1);
    
    % modifiy edit field value
    set_non_nan_value(h.phaseUnwrap.edit.excludeMask, 'String', val);
    
    % popup manu for thresholding method
    str_pattern	= '.unwrap.excludeMethod';
    val      	= get_string_as_string(config_txt, str_pattern);
    if ~isnan(val)
        % matching algorithm name
        for k = 1:length(methodExcludedName)
            if strcmpi(val,methodExcludedName{k})
                val = k;
                break
            end
        end
    
        % change method popup manu
        set_non_nan_value(h.phaseUnwrap.popup.excludeMethod,'Value',val);
    end
end

% save unwrap
str_pattern	= '.unwrap.isSaveUnwrappedEcho';
val      	= get_num_as_string(config_txt, str_pattern,'=',';');
set_non_nan_value(h.phaseUnwrap.checkbox.saveEchoPhase, 'Value', str2double(val));

end

%% Panel: BFR
if isModifyBFR
    
% % logical B1 polyfit
% str_pattern	= '.bfr.refine';
% val      	= get_num_as_string(config_txt, str_pattern,'=',';');
% set_non_nan_value(h.bkgRemoval.checkbox.refine, 'Value', str2double(val));
% B1 polyfit method
% popup manu for fit method
str_pattern	= '.bfr.refine_method';
val      	= get_string_as_string(config_txt, str_pattern);
if ~isnan(val)
    % matching algorithm name
    for k = 1:length(methodRefineName)
        if strcmpi(val,methodRefineName{k})
            val = k;
            break
        end
    end

    % change method popup manu
    set_non_nan_value(h.bkgRemoval.popup.refine,'Value',val);
    
    % trigger popup callback 
    feval(h.bkgRemoval.popup.refine.Callback{1},h.bkgRemoval.popup.refine,[],h)

end
% B1 polyfit order
str_pattern	= '.bfr.refine_order';
val      	= get_num_as_string(config_txt, str_pattern,'=',';');
set_non_nan_value(h.bkgRemoval.edit.order, 'String', val);

% Erosion radius
str_pattern	= '.bfr.erode_radius';
val      	= get_num_as_string(config_txt, str_pattern,'=',';');
set_non_nan_value(h.bkgRemoval.edit.imerode, 'String', val);

% method
str_pattern	= '.bfr.method';
BFR_method	= get_string_as_string(config_txt, str_pattern);
if ~isnan(BFR_method)
    
    for k = 1:length(methodBFRName)
        if strcmpi(BFR_method,methodBFRName{k})
            feval(config_BFR_function{k},h,'get',config_txt);
            BFR_method = k;
        end
    end
    
    % change method popup manu
    set_non_nan_value(h.bkgRemoval.popup.bkgRemoval,'Value',BFR_method);
    % trigger popup callback to switch method panel
    feval(h.bkgRemoval.popup.bkgRemoval.Callback{1},h.bkgRemoval.popup.bkgRemoval,[],h)
end

end
%% Panel: QSM
if isModifyQSM
    
% reference tissue
str_pattern	= '.qsm.reference_tissue';
val      	= get_string_as_string(config_txt, str_pattern);
switch lower(val)
    case 'none'
        set_non_nan_value(h.qsm.popup.tissue,'Value',1)
    case 'brain mask'
        set_non_nan_value(h.qsm.popup.tissue,'Value',2)
    case 'csf'
        set_non_nan_value(h.qsm.popup.tissue,'Value',3)
end

% method
str_pattern	= '.qsm.method';
QSM_method	= get_string_as_string(config_txt, str_pattern);
if ~isnan(QSM_method)
    
    for k = 1:length(methodQSMName)
        if strcmpi(QSM_method,methodQSMName{k})
            feval(config_QSM_function{k},h,'get',config_txt);
            QSM_method = k;
        end
    end
    
    % change method popup manu
    set_non_nan_value(h.qsm.popup.qsm,'Value',QSM_method);
    % trigger popup callback to switch method panel
    feval(h.qsm.popup.qsm.Callback{1},h.qsm.popup.qsm,[],h)
end

end

end
