%% sepia_read_popup_value(config_txt, str_pattern, action_handle, popup_list)
%
% Input
% --------------
% config_filename : filename of the sepia pipeline configureation file
% h               : master handle of SEPIA
%
% Description: Read pipeline configuration in the file back to GUI
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 6 March 2020 (v0.8.1)
% Date modified:12 June 2021 (v1.0)
% Date modified:5 July 2025 (v1.3)
%
%
function set_config_Callback(config_filename,h)

sepia_universal_variables;

config_txt = fileread(config_filename);

%% to determine what is going to be modified
isModifyGeneral = false;
isModifyUnwrap  = false;
isModifyBFR     = false;
isModifyQSM     = false;
isModifyIO      = false;

str_pattern	= 'input';
idx         = regexp(config_txt,str_pattern,'once');
if ~isempty(idx)
    isModifyIO = true;
end
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
if isModifyIO
input_action_handle  = h.dataIO.edit.input;
input1_action_handle = h.dataIO.edit.inputData1;
input2_action_handle = h.dataIO.edit.inputData2;
input3_action_handle = h.dataIO.edit.inputData3;
header_action_handle = h.dataIO.edit.inputHeader;
mask_action_handle   = h.dataIO.edit.maskdir;
output_action_handle = h.dataIO.edit.output;

isAllInputEmpty = (isempty(input_action_handle.String) .* isempty(input1_action_handle.String) .* ...
                  isempty(input2_action_handle.String).* isempty(input3_action_handle.String) .* ...
                  isempty(header_action_handle.String))>0;

isModifyOutput = true;

% check which input method was used in the config file
ind = regexp(config_txt, regexptranslate('wildcard','input(*).name'), 'once');
if ~isempty(ind)
    isFileIO = true;
else
    isFileIO = false;
end

% Directory input method
action_handle = input_action_handle;
if isAllInputEmpty && ~isFileIO     % only modify the field if all IO fields are empty
    str_pattern     = 'input';
    str_end_idx     = regexp(config_txt,str_pattern,'end');
    indicator_idx   = regexp(config_txt,'''');
    str = config_txt(indicator_idx(find(indicator_idx > str_end_idx(1), 1 ))+1:indicator_idx(find(indicator_idx > str_end_idx(1), 1 )+1)-1);

    if prod(~isnan(str)) && isfolder(str)
        set(action_handle,'String',str);
    else
        isModifyOutput = false;
    end
end

% File input method
% input 1
action_handle = input1_action_handle;
if and(isempty(input_action_handle.String),isempty(get(action_handle,'String')))
    str_pattern     = 'input\(1).name';
    str = get_string_as_string(config_txt, str_pattern);
    if prod(~isnan(str)) && isfile(str)
        set(action_handle,'String',str);
    else
        isModifyOutput = false;
    end
end

% input 2
action_handle = input2_action_handle;
if and(isempty(input_action_handle.String),isempty(get(action_handle,'String')))
    str_pattern     = 'input\(2).name';
    str = get_string_as_string(config_txt, str_pattern);
    if prod(~isnan(str)) && isfile(str)
        set(action_handle,'String',str);
    else
        isModifyOutput = false;
    end
end

% input 3
action_handle = input3_action_handle;
if and(isempty(input_action_handle.String),isempty(get(action_handle,'String')))
    str_pattern     = 'input\(3).name';
    str = get_string_as_string(config_txt, str_pattern);
    if prod(~isnan(str)) && isfile(str)
        set(action_handle,'String',str);
    else
        isModifyOutput = false;
    end
end

% header input
action_handle = header_action_handle;
if and(isempty(input_action_handle.String),isempty(get(action_handle,'String')))
    str_pattern     = 'input\(4).name';
    str = get_string_as_string(config_txt, str_pattern);
    if prod(~isnan(str)) && isfile(str)
        set(action_handle,'String',str);
    else
        isModifyOutput = false;
    end
end

% mask 
action_handle = mask_action_handle;
if isempty(get(action_handle,'String'))
    str_pattern     = 'mask_filename = ';
    str = get_string_as_string(config_txt, str_pattern);
    if prod(~isnan(str))
        set(action_handle,'String',str);
    end
end

% output 
action_handle = output_action_handle;
if isempty(get(action_handle,'String'))
    str_pattern     = 'output_basename = ';
    str = get_string_as_string(config_txt, str_pattern);
    if prod(~isnan(str)) && isModifyOutput
        set(action_handle,'String',str);
    end
end


end

if isModifyGeneral
    
% invert phase data
str_pattern     = '.general.isInvert';
action_handle   = h.dataIO.checkbox.invertPhase;
sepia_read_checkbox_value(config_txt, str_pattern, action_handle);

% BET
str_pattern     = '.general.isBET';
action_handle   = h.dataIO.checkbox.brainExtraction;
val = sepia_read_checkbox_value(config_txt, str_pattern, action_handle);
% trigger checkout callback
feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);

% isBet is true then change the BET parameters
if str2double(val)

    % 20250705 v1.3
    str_pattern     = '.general.brain_extraction_method';
    action_handle   = h.dataIO.popup.brainExtraction;
    sepia_read_popup_value(config_txt, str_pattern, action_handle, skullstrippingMethod);
    % trigger checkout callback
    feval(h.dataIO.popup.brainExtraction.Callback{1},h.dataIO.popup.brainExtraction,[],h);

    % -f
    str_pattern     = '.general.fractional_threshold';
    action_handle   = h.dataIO.edit.fractionalThres;
    sepia_read_edit_string(config_txt, str_pattern, action_handle);
    
    % -g
    str_pattern     = '.general.gradient_threshold';
    action_handle   = h.dataIO.edit.gradientThres;
    sepia_read_edit_string(config_txt, str_pattern, action_handle);
    
end

% refine brain mask
str_pattern     = '.general.isRefineBrainMask';
action_handle   = h.dataIO.checkbox.refineBrainMask;
sepia_read_checkbox_value(config_txt, str_pattern, action_handle);

% 20250705 v1.3 popup manu for denoise
str_pattern     = '.general.isDenoise';
action_handle   = h.dataIO.checkbox.denoise;
val = sepia_read_checkbox_value(config_txt, str_pattern, action_handle);
% trigger callback 
feval(action_handle.Callback{1},action_handle,[],{h.dataIO.edit.denoise,h.dataIO.slider.denoise},1); 
if str2double(val)
% modifiy edit field value
str_pattern     = '.general.denoiseKernel';
val             = get_num_as_string(config_txt, str_pattern, '=', ';');
set_non_nan_value(h.dataIO.edit.denoise, 'String', val);
end
% 20250705 v1.3 popup manu for denoise
str_pattern     = '.general.isUpsample';
action_handle   = h.dataIO.checkbox.upsample;
val = sepia_read_checkbox_value(config_txt, str_pattern, action_handle);
% trigger callback 
feval(action_handle.Callback{1},action_handle,[],{h.dataIO.edit.upsample,h.dataIO.slider.upsample},1); 
if str2double(val)
% modifiy edit field value
str_pattern     = '.general.target_resolution';
val             = get_num_as_string(config_txt, str_pattern, '=', ';');
set_non_nan_value(h.dataIO.edit.upsample, 'String', val);
end

end

%% Panel: Total field
if isModifyUnwrap
    
% echo combination method
str_pattern      = '.unwrap.echoCombMethod';
action_handle    = h.phaseUnwrap.popup.phaseCombMethod;
popup_list       = methodEchoCombineName;
config_func_list = config_EchoCombine_function;
read_method_popup(config_txt, str_pattern, action_handle, popup_list, config_func_list, h)

end

%% Panel: BFR
if isModifyBFR
    
% B1 polyfit method
% popup manu for fit method
str_pattern     = '.bfr.refine_method';
action_handle   = h.bkgRemoval.popup.refine;
sepia_read_popup_value(config_txt, str_pattern, action_handle, methodRefineName);
% trigger popup callback 
feval(h.bkgRemoval.popup.refine.Callback{1},h.bkgRemoval.popup.refine,[],h)

% B1 polyfit order
str_pattern     = '.bfr.refine_order';
action_handle   = h.bkgRemoval.edit.order;
sepia_read_edit_string(config_txt, str_pattern, action_handle);

% Erosion radius
str_pattern     = '.bfr.erode_radius';
action_handle   = h.bkgRemoval.edit.imerode;
sepia_read_edit_string(config_txt, str_pattern, action_handle);

% Erosion radius before BFR
try % backward compatible
str_pattern     = '.bfr.erode_before_radius';
action_handle   = h.bkgRemoval.edit.imerodebefore;
sepia_read_edit_string(config_txt, str_pattern, action_handle);
catch
end

% background field removal method
str_pattern      = '.bfr.method';
action_handle    = h.bkgRemoval.popup.bkgRemoval;
popup_list       = methodBFRName;
config_func_list = config_BFR_function;
read_method_popup(config_txt, str_pattern, action_handle, popup_list, config_func_list, h)

end
%% Panel: QSM
if isModifyQSM
    
% reference tissue
str_pattern     = '.qsm.reference_tissue';
action_handle   = h.qsm.popup.tissue;
sepia_read_popup_value(config_txt, str_pattern, action_handle, tissueName);

% QSM dipole inversion method
str_pattern      = '.qsm.method';
action_handle    = h.qsm.popup.qsm;
popup_list       = methodQSMName;
config_func_list = config_QSM_function;
read_method_popup(config_txt, str_pattern, action_handle, popup_list, config_func_list, h)

end

end

%% read the popup method that triggers the switch of method panel
% config_txt    : variable contains config text
% str_pattern   : string pattern to be printed after algorParam parameter
% action_handle : handle of the GUI popup
% popup_list    : popup list, in cell (methods in this case)
% config_func_list: the correspinding get_set_XXX file of the methods
% h             : master handle
function read_method_popup(config_txt, str_pattern, action_handle, popup_list, config_func_list, h)

% get option as string
method = get_string_as_string(config_txt, str_pattern);

if ~isnan(method)
    % matching popup list name
    for j = 1:length(popup_list)
        if strcmpi(method,popup_list{j})
            method = j;
            feval(config_func_list{j},h,'get',config_txt);
            break
        end
    end

    % change popup manu
    set_non_nan_value(action_handle,'Value',method);
    % trigger popup callback to switch method panel
    feval(action_handle.Callback{1},action_handle,[],h)
end

end
