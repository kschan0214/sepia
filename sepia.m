%% sepia
%
% Description: This is a GUI of sepia, which is a pipeline control tool
% for standard QSM processing. It supports the following processing steps:
% (1) phase unwrapping
% (2) background field removal
% (3) QSM
% All these processing steps can be performed in one go or worked as
% standalones.
%
% After perfoming all the selected processing steps sucessfully, a log file
% will be generated. This log file contains the text of how a function is
% called to run the process(es). You can also use the content of the log
% file to perform your data processing without calling this GUI  
% (which in some cases good for batch job).
%
% To run this GUI, simply add the directory containing this file in your
% Matlab's PATH. This GUI will automatically add the related files into
% your Matlab's PATH during processing.
%
% Considering different toolboxes may use their own function for the same
% processing which has the same name, the specific method will only be
% loaded just before being used as Matlab does not support function
% overloading. Nevertheless the GUI will automatically do this for you.
%
% Most of the processing methods are provided from other toolboxes,
% inclduing MEDI, STI Suite and FANSI. Please cite the correspoding
% paper(s) when you use this GUI.
%
% **********************************
% * Hierarchy of the GUI structure *
% **********************************
% |---fig
%   |---Tabs
%     |---StepsPanel
%       |---Method(also Edit,Text,etc.)
%         |---Edit, text, Popup, etc.
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 September 2017
% Date modified: 1 June 2018
% Date modified: 2 March 2020 (v0.8.0)
% Date modified: 27 Jan 2020 (v0.8.1)
% Date modified: 12 June 2021 (v1.0)
% Date modified: 3 August 2022 (v1.1)
% Date modified: 3 April 2023 (v1.2.2.4)
% Date modified: 9 October 2023 (v1.2.2.5)
% Date modified: 7 July 2025 (v1.3)
%
function sepia 

% clear previous handles
clear global h 

% make sure nothing is logged at the moment
diary off

% add path and check toolboxes availability
sepia_addpath('',1);

global h 

% SEPIA version
sepia_universal_variables;
% SEPIA_version = 'v0.8.1';

%% create basic GUI
% set GUI window size
screenSize  = get(0,'ScreenSize');
posLeft     = round(screenSize(3)/5);
posBottom   = round(screenSize(4)/8);
guiSizeHori = round(screenSize(3)/1.5);
guiSizeVert = round(screenSize(4)*3/4);
if guiSizeHori < 1000
    guiSizeHori = 1000;
end
if guiSizeVert < 700
    guiSizeVert = 700;
end

% create GUI figure
h.fig=figure('Units','pixels','position',[posLeft posBottom guiSizeHori guiSizeVert],...
    'MenuBar','None','Toolbar','None','Name',['SEPIA GUI (' SEPIA_version ')'],'NumberTitle','off');

%% construct panels for each tab 
% create Tabs for GUI
h.TabGroup          = uitabgroup(h.fig,'position',[.01 .01 0.98 0.98]);
h.Tabs.Sepia        = uitab(h.TabGroup,'Title','SEPIA');
h.Tabs.phaseUnwrap  = uitab(h.TabGroup,'Title','Phase unwrapping');
h.Tabs.bkgRemoval   = uitab(h.TabGroup,'Title','Background field removal');
h.Tabs.qsm          = uitab(h.TabGroup,'Title','QSM');
h.Tabs.swismwi      = uitab(h.TabGroup,'Title','SWI/SMWI');
h.Tabs.r2star       = uitab(h.TabGroup,'Title','R2* mapping');
h.Tabs.analysis   	= uitab(h.TabGroup,'Title','Analysis');
h.Tabs.utility      = uitab(h.TabGroup,'Title','Utility');

%% GUI with QSM one-stop station tab
% these panels will be switching position from tab to tab
parent_curr = h.Tabs.Sepia;
% I/O
h = sepia_handle_panel_dataIO(parent_curr,      h,[0.01 0.8]);
% phase unwrap
h = sepia_handle_panel_phaseUnwrap(parent_curr,	h,[0.01 0.59]);
% background field
h = sepia_handle_panel_bkgRemoval(parent_curr,	h,[0.01 0.33]);
% QSM
h = sepia_handle_panel_qsm(parent_curr,       	h,[0.01 0.07]);

%% SWI/SMWI tab
parent_curr = h.Tabs.swismwi;
% I/O
h = sepia_handle_panel_swismwi_dataIO(parent_curr,	h,[0.01 0.8]);
% Method
h = sepia_handle_panel_swismwi(parent_curr,       	h,[0.01 0.44]);

%% R2* mapping tab
parent_curr = h.Tabs.r2star;
% R2*
h = sepia_handle_panel_r2s(parent_curr,	h,[0.01 0.55]);

%% analysis tab
h = sepia_handle_panel_Analysis(h.Tabs.analysis,	h,[0.01 0.01]);

%% utility tab
h = sepia_handle_panel_Utility(h.Tabs.utility,      h,[0.01 0.39]);

%% extra content
% Start button
h.pushbutton_start = uicontrol('Parent',h.Tabs.Sepia,...
    'Style','pushbutton',...
    'String','Start',...
    'units','normalized','Position',[0.85 0.01 0.1 0.05],...
    'backgroundcolor',get(h.fig,'color'));

% load config button
h.pushbutton_loadConfig = uicontrol('Parent',h.Tabs.Sepia,...
    'Style','pushbutton',...
    'String','Load config',...
    'units','normalized','Position',[0.01 0.01 0.1 0.05],...
    'backgroundcolor',get(h.fig,'color'));

%% Set Callback functions
set(h.TabGroup,                 'SelectionChangedFcn', {@SwitchTab_Callback})
set(h.pushbutton_start,         'Callback',            {@PushbuttonStart_Callback});
set(h.pushbutton_loadConfig,    'Callback',            {@PushbuttonLoadConfig_Callback});

end

%% Callback functions

%% switching tabs
function SwitchTab_Callback(source,eventdata)
% switch parent handle of StepsPanel based on current tab

global h tooltip fieldString

% global uicontrol for most of the tabs
universial_handle = {h.StepsPanel.dataIO,...
                     h.pushbutton_start,...         % Start pushbutton
                     h.pushbutton_loadConfig,...    % load config pushbutton
                     };               % GPU checkbox

% Tab specific strings and tooltips
tooltip.input_dir{1} = 'Directory contains phase (*ph*.nii*), magnitude (*mag*.nii) & header (*header*.mat) (& mask, *mask*nii*) files';
tooltip.input_dir{2} = 'Directory contains the total field map (*fieldmap*.nii*) and SEPIA header (*header*.mat)';
tooltip.input_dir{3} = 'Directory contains the local field map (*localfield*.nii*) (depending on QSM algorithm, additional file(s) may also be needed, e.g. *mag*.nii* and *weights*.nii*)';

fieldString.inputData1{1}= 'or Phase:';
fieldString.inputData1{2}= 'or Fieldmap:';
fieldString.inputData1{3}= 'or Local field:';

fieldString.inputData3{1}= '    Weights:';
fieldString.inputData3{2}= '    Noise SD:';

% change universal elements' parent except 'SWI/SMWI' tab and 'Utility' tab
if ~strcmpi(eventdata.NewValue.Title,'SWI/SMWI') && ~strcmpi(eventdata.NewValue.Title,'Utility') && ~strcmpi(eventdata.NewValue.Title,'Analysis')
    for k = 1:length(universial_handle)
        set(universial_handle{k}, 'Parent', source.SelectedTab);
    end
end

% Specify tab-specific content
switch eventdata.NewValue.Title
    
    % QSM one-stop station tab
    case 'SEPIA'
        switch_tab_to_SEPIA;

    % Phase unwrapping tab
    case 'Phase unwrapping'
        switch_tab_to_phase_unwrapping;
        
    % background field removal tab    
    case 'Background field removal'
        switch_tab_to_BRF;
        
    % qsm tab    
    case 'QSM'
        switch_tab_to_QSM;

    case 'R2* mapping'
        switch_tab_to_R2s;
        
    case 'Analysis'
        switch_tab_to_Analysis;

end

end

%% Start button
function PushbuttonStart_Callback(source,eventdata)
% objective of Sepia GUI is to create a .m script to execute a command-based function

global h

sepia_universal_variables;

%%%%%%%%%%%% Step 1: preparation %%%%%%%%%%%%
% Disable the pushbutton to prevent double clicks
set(source,'Enable','off');

% get the current tab name of the standalone
tab = h.StepsPanel.dataIO.Parent.Title;

% get I/O GUI input
% option 1: input is a directory
input           = get(h.dataIO.edit.input,          'String');
% option 2: input are NIfTI files
if isempty(input)
    input(1).name = get(h.dataIO.edit.inputData1, 	'String');
    input(2).name = get(h.dataIO.edit.inputData2,  	'String');
    input(3).name = get(h.dataIO.edit.inputData3,  	'String');
    input(4).name = get(h.dataIO.edit.inputHeader, 	'String');
end
outputBasename  = get(h.dataIO.edit.output,       	'String');
maskFullName    = get(h.dataIO.edit.maskdir,       	'String');

% get and create output directory
output_index    = strfind(outputBasename, filesep);
outputDir       = outputBasename(1:output_index(end));
% if the output directory does not exist then create the directory
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end


% use current time as unique identifier
identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');
% create a new m file
configFilename = fullfile(outputDir, ['sepia_config' identifier '.m']);
% configFilename = [outputDir filesep 'sepia_config.m'];
if exist(configFilename,'file') == 2
    % get new time index
    identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');
    configFilename = fullfile(outputDir, ['sepia_config' identifier '.m']);

end
fid = fopen(configFilename,'w');

%%%%%%%%%%%% Step 2: Write config file %%%%%%%%%%%% 
% specify config file version
fprintf(fid,'%% This file is generated by SEPIA version: %s\n',SEPIA_version);
% general path
fprintf(fid,'%% add general Path\n');
fprintf(fid,'sepia_addpath;\n\n');

fprintf(fid,'%% Input/Output filenames\n');
fprintf(fid,'input = struct();\n');
% input data
if isstruct(input)
    fprintf(fid,'input(1).name = ''%s'' ;\n',input(1).name);
    fprintf(fid,'input(2).name = ''%s'' ;\n',input(2).name);
    fprintf(fid,'input(3).name = ''%s'' ;\n',input(3).name);
    fprintf(fid,'input(4).name = ''%s'' ;\n',input(4).name);
else
    fprintf(fid,'input = ''%s'' ;\n',input);
end

% output
fprintf(fid,'output_basename = ''%s'' ;\n',outputBasename);

% mask
fprintf(fid,'mask_filename = [''%s''] ;\n\n',maskFullName);

% general algorithm parameters
fprintf(fid,'%% General algorithm parameters\n');
fprintf(fid,'algorParam = struct();\n');
% BET
isbet = sepia_print_checkbox_value(fid,'.general.isBET',h.dataIO.checkbox.brainExtraction);
sepia_print_popup_as_string(fid, '.general.brain_extraction_method',h.dataIO.popup.brainExtraction);
if isbet
    switch h.dataIO.popup.brainExtraction.String{h.dataIO.popup.brainExtraction.Value,1} 
        case 'FSL bet (MEDI)'
            % for backward compatability
            sepia_print_edit_as_string(fid,'.general.fractional_threshold',h.dataIO.edit.fractionalThres);
            sepia_print_edit_as_string(fid,'.general.gradient_threshold',h.dataIO.edit.gradientThres);

    end
    
end
sepia_print_checkbox_value(fid,'.general.isInvert',h.dataIO.checkbox.invertPhase);

% refine brain mask
sepia_print_checkbox_value(fid,'.general.isRefineBrainMask',h.dataIO.checkbox.refineBrainMask);

% denoise
isDenoise = sepia_print_checkbox_value(fid,'.general.isDenoise',h.dataIO.checkbox.denoise);
if isDenoise
    sepia_print_edit_as_string(fid,'.general.denoiseKernel',h.dataIO.edit.denoise);
end
% upsample
isUpsample = sepia_print_checkbox_value(fid,'.general.isUpsample',h.dataIO.checkbox.upsample);
if isUpsample
    sepia_print_edit_as_string(fid,'.general.target_resolution',h.dataIO.edit.upsample);
end

% phase unwrap algorithm parameters
if strcmpi(tab,'SEPIA') || strcmpi(tab,'Phase unwrapping')
    fprintf(fid,'%% Total field recovery algorithm parameters\n');
    
    % echo phase combine method
    % set parameters for selected method
    print_method_popup_and_eval(fid, '.unwrap.echoCombMethod', h.phaseUnwrap.popup.phaseCombMethod, methodEchoCombineName, config_EchoCombine_function, h);

end
                
% background field removal algorithm parameters
if strcmpi(tab,'SEPIA') || strcmpi(tab,'Background field removal')
    fprintf(fid,'%% Background field removal algorithm parameters\n');
    % polyfit method
    sepia_print_popup_as_string(fid,'.bfr.refine_method',h.bkgRemoval.popup.refine);
    % polyfit order
    sepia_print_edit_as_string(fid,'.bfr.refine_order',h.bkgRemoval.edit.order);
    % Erode local field
    sepia_print_edit_as_string(fid,'.bfr.erode_radius',h.bkgRemoval.edit.imerode);
    % Erode mask before BFR
    sepia_print_edit_as_string(fid,'.bfr.erode_before_radius',h.bkgRemoval.edit.imerodebefore);
    
    % set parameters for selected method
    print_method_popup_and_eval(fid, '.bfr.method', h.bkgRemoval.popup.bkgRemoval, methodBFRName, config_BFR_function, h);

end

% QSM algorithm parameters
if strcmpi(tab,'SEPIA') || strcmpi(tab,'QSM')
    fprintf(fid,'%% QSM algorithm parameters\n');
    
    % reference tissue
    sepia_print_popup_as_string(fid,'.qsm.reference_tissue',h.qsm.popup.tissue);
    
    % set parameters for selected method
    print_method_popup_and_eval(fid, '.qsm.method', h.qsm.popup.qsm, methodQSMName, config_QSM_function, h);

    % HEIDI
    isHEIDI = sepia_print_checkbox_value(fid,'.qsm.isHEIDI',h.qsm.checkbox.isHeidi);
    if isHEIDI
        if ~strcmp(h.qsm.popup.qsm.String{h.qsm.popup.qsm.Value},'LSQR+HEIDI')
            print_HEIDI_popup_and_eval(fid, methodQSMName, config_QSM_function, h);
        else
            warning('The selected dipole inversion method already includes HEIDI. No extra HEIDI processing will be done.');
        end
    end

end

% R2* algorithm parameters
if strcmpi(tab,'R2* mapping')
    fprintf(fid,'%% R2* algorithm parameters\n');
    
    % set parameters for selected method
    print_method_popup_and_eval(fid, '.r2s.method', h.r2s.popup.r2sMethod, methodR2sName, config_R2s_function, h);

end

% Determine application based on Tab
fprintf(fid,'\nsepiaIO(input,output_basename,mask_filename,algorParam);\n');

fclose(fid);

try
    % run process
    run(configFilename);
    
%     % turn off the log
%     diary off
    
    % re-enable the pushbutton
    set(source,'Enable','on');

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    
    % rethrow the error message to command window
    rethrow(ME);
    
end


end

function PushbuttonLoadConfig_Callback(source,eventdata)

global h tooltip fieldString

% only read m file
[config_filename,pathDir] = uigetfile({'*.m'},'Select a SEPIA config file');

% if file is specified then read it
if exist(fullfile(pathDir,config_filename),'file')
    set_config_Callback(fullfile(pathDir,config_filename),h);
end

tab = h.StepsPanel.dataIO.Parent.Title;

% Tab specific strings and tooltips
tooltip.input_dir{1} = 'Directory contains phase (*ph*.nii*), magnitude (*mag*.nii) & header (*header*.mat) (& mask, *mask*nii*) files';
tooltip.input_dir{2} = 'Directory contains the total field map (*fieldmap*.nii*) and SEPIA header (*header*.mat)';
tooltip.input_dir{3} = 'Directory contains the local field map (*localfield*.nii*) (depending on QSM algorithm, additional file(s) may also be needed, e.g. *mag*.nii* and *weights*.nii*)';

fieldString.inputData1{1}= 'or Phase:';
fieldString.inputData1{2}= 'or Fieldmap:';
fieldString.inputData1{3}= 'or Local field:';

fieldString.inputData3{1}= '    Weights:';
fieldString.inputData3{2}= '    Noise SD:';

switch tab
    
    % QSM one-stop station tab
    case 'SEPIA'
        switch_tab_to_SEPIA;
        
    % Phase unwrapping tab
    case 'Phase unwrapping'
        switch_tab_to_phase_unwrapping;
        
    % background field removal tab    
    case 'Background field removal'
        switch_tab_to_BRF;
    % qsm tab    
    case 'QSM'
        switch_tab_to_QSM

    case 'R2* mapping'
        switch_tab_to_R2s
        
end

end

function print_method_popup_and_eval(fid, str_pattern, action_handle, popup_list, config_list, h)

% set parameters for selected method
method = sepia_print_popup_as_string(fid,str_pattern,action_handle);
for k = 1:length(popup_list)
    if strcmpi(method,popup_list{k})
        feval(config_list{k},h,'set',fid);
    end
end
    
end

function print_HEIDI_popup_and_eval(fid, popup_list, config_list, h)

% set parameters for selected method
for k = 1:length(popup_list)
    if strcmpi('LSQR+HEIDI',popup_list{k})
        feval(config_list{k},h,'set',fid);
    end
end
    
end

function switch_tab_to_SEPIA
    global h tooltip fieldString
    % I/O
    % Change essential files if input is a directory
    set(h.dataIO.text.input,                'Tooltip',tooltip.input_dir{1});
    % BET is supported with this tab
    set(h.dataIO.checkbox.brainExtraction,  'Enable','on');
    set(h.dataIO.popup.brainExtraction,     'Enable','on');
        % trigger followup callback to switch method panel
        feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
    % phase invert is supported with this tab
    set(h.dataIO.checkbox.invertPhase,      'Enable','on');
    % input data 1
    set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{1});
    set(h.dataIO.edit.inputData1,           'Enable','on');
    set(h.dataIO.button.inputData1,         'Enable','on');
    % input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
    set(h.dataIO.edit.inputData2,           'Enable','on');
    set(h.dataIO.button.inputData2,         'Enable','on');
    % input data 3
    set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{1});
    set(h.dataIO.edit.inputData3,           'Enable','on');
    set(h.dataIO.button.inputData3,         'Enable','on');
    % refine brain mask is supported with this tab
    set(h.dataIO.checkbox.refineBrainMask,  'Enable','on');
    % denoise and upsample
    set(h.dataIO.checkbox.denoise,          'Enable','on');
    set(h.dataIO.checkbox.upsample,         'Enable','on');

    % phase unwrap
    set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.Sepia,'Position',[0.01 0.59 0.95 0.2]);
    % background field
    set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.Sepia,'Position',[0.01 0.33 0.95 0.25]);
    % QSM
    set(h.StepsPanel.qsm,           'Parent',h.Tabs.Sepia,'Position',[0.01 0.07 0.95 0.25]);
end

function switch_tab_to_phase_unwrapping
global h tooltip fieldString

% I/O
% This tab supports both DICOM and NIfTI files
set(h.dataIO.text.input,                'Tooltip',tooltip.input_dir{1});
% BET is supported with this tab
set(h.dataIO.checkbox.brainExtraction,  'Enable','on');
set(h.dataIO.popup.brainExtraction,     'Enable','on');
    % trigger followup callback to switch method panel
    feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
% phase invert is supported with this tab
set(h.dataIO.checkbox.invertPhase,      'Enable','on');
% input data 1
set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{1});
set(h.dataIO.edit.inputData1,           'Enable','on');
set(h.dataIO.button.inputData1,         'Enable','on');
% input data 2
%           set(h.dataIO.text.inputData2,'String','Magn. data:');
set(h.dataIO.edit.inputData2,           'Enable','on');
set(h.dataIO.button.inputData2,         'Enable','on');
% input data 3
set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{1});
set(h.dataIO.edit.inputData3,           'Enable','off','String',[]);
set(h.dataIO.button.inputData3,         'Enable','off');
% refine brain mask is supported with this tab
set(h.dataIO.checkbox.refineBrainMask,  'Enable','on');
% denoise and upsample
set(h.dataIO.checkbox.denoise,          'Enable','on');
set(h.dataIO.checkbox.upsample,         'Enable','on');

% phase unwrap
set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.phaseUnwrap,'Position',[0.01 0.59 0.95 0.2]);
end

function switch_tab_to_BRF
global h tooltip fieldString

% I/O
% This tab supports only NIfTI files
set(h.dataIO.text.input,                'Tooltip',tooltip.input_dir{2});
% no BET support with this tab
set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
set(h.dataIO.popup.brainExtraction,     'Enable','on');
    % trigger followup callback to switch method panel
    feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
set(h.dataIO.edit.maskdir,              'Enable','on');
set(h.dataIO.button.maskdir,            'Enable','on');
% phase invert is not supported with this tab
set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
% input data 1
set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{2});
set(h.dataIO.edit.inputData1,           'Enable','on');
set(h.dataIO.button.inputData1,         'Enable','on');
% input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
set(h.dataIO.edit.inputData2,           'Enable','off','String',[]);
set(h.dataIO.button.inputData2,         'Enable','off');
% input data 3
set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{2});
set(h.dataIO.edit.inputData3,           'Enable','on');
set(h.dataIO.button.inputData3,         'Enable','on');
% no refine brain mask with this tab
set(h.dataIO.checkbox.refineBrainMask,  'Enable','off','Value',0);
% denoise and upsample
set(h.dataIO.checkbox.denoise,          'Enable','off','Value',0);
set(h.dataIO.checkbox.upsample,         'Enable','off','Value',0);
    % trigger followup callback to switch method panel
    feval(h.dataIO.checkbox.denoise.Callback{1},h.dataIO.checkbox.denoise,[],{h.dataIO.edit.denoise,h.dataIO.slider.denoise},1);
    feval(h.dataIO.checkbox.upsample.Callback{1},h.dataIO.checkbox.upsample,[],{h.dataIO.edit.upsample,h.dataIO.slider.upsample},1);

% background field
set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.bkgRemoval,'Position',[0.01 0.54 0.95 0.25]);

end

function switch_tab_to_QSM
global h tooltip fieldString
% I/O
% This tab supports only NIfTI files
set(h.dataIO.text.input,                'Tooltip',tooltip.input_dir{3});
% no BET support with this tab
set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
set(h.dataIO.popup.brainExtraction,     'Enable','off');
    % trigger followup callback to switch method panel
    feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
set(h.dataIO.edit.maskdir,              'Enable','on');
set(h.dataIO.button.maskdir,            'Enable','on');
% phase invert is not supported with this tab
set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
% input data 1
set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{3});
set(h.dataIO.edit.inputData1,           'Enable','on');
set(h.dataIO.button.inputData1,         'Enable','on');
% input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
set(h.dataIO.edit.inputData2,           'Enable','on');
set(h.dataIO.button.inputData2,         'Enable','on');
% input data 3
set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{1});
set(h.dataIO.edit.inputData3,           'Enable','on');
set(h.dataIO.button.inputData3,         'Enable','on');
% no refine brain mask with this tab
set(h.dataIO.checkbox.refineBrainMask,  'Enable','off','Value',0);
% denoise and upsample
set(h.dataIO.checkbox.denoise,          'Enable','off','Value',0);
set(h.dataIO.checkbox.upsample,         'Enable','off','Value',0);
    % trigger followup callback to switch method panel
    feval(h.dataIO.checkbox.denoise.Callback{1},h.dataIO.checkbox.denoise,[],{h.dataIO.edit.denoise,h.dataIO.slider.denoise},1);
    feval(h.dataIO.checkbox.upsample.Callback{1},h.dataIO.checkbox.upsample,[],{h.dataIO.edit.upsample,h.dataIO.slider.upsample},1);
% QSM
set(h.StepsPanel.qsm,           'Parent',h.Tabs.qsm,'Position',[0.01 0.54 0.95 0.25]);
end

function switch_tab_to_R2s
global h tooltip fieldString
% BET is not supported with this tab
set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
    % trigger followup callback to switch method panel
    feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
set(h.dataIO.edit.maskdir,              'Enable','on');
set(h.dataIO.button.maskdir,            'Enable','on');
% phase invert is not supported with this tab
set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
% no refine brain mask with this tab
set(h.dataIO.checkbox.refineBrainMask,  'Enable','off','Value',0);
% input data 1
set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{1});
set(h.dataIO.edit.inputData1,           'Enable','off','String',[]);
set(h.dataIO.button.inputData1,         'Enable','off');
% input data 3
set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{1});
set(h.dataIO.edit.inputData3,           'Enable','off','String',[]);
set(h.dataIO.button.inputData3,         'Enable','off');
% denoise and upsample
set(h.dataIO.checkbox.denoise,          'Enable','off','Value',0);
set(h.dataIO.checkbox.upsample,         'Enable','off','Value',0);
    % trigger followup callback to switch method panel
    feval(h.dataIO.checkbox.denoise.Callback{1},h.dataIO.checkbox.denoise,[],{h.dataIO.edit.denoise,h.dataIO.slider.denoise},1);
    feval(h.dataIO.checkbox.upsample.Callback{1},h.dataIO.checkbox.upsample,[],{h.dataIO.edit.upsample,h.dataIO.slider.upsample},1);
end

function switch_tab_to_Analysis
% test if the directory exist
try 
    SpecifyToolboxesDirectory;
    if exist('ANTS_HOME', 'var')
        if ~exist(ANTS_HOME,'dir')
            warndlg('Missing ANTs library. All functions in the Analysis Tab cannot be used.')
        end
    else
            warndlg('Missing ANTs library. All functions in the Analysis Tab cannot be used.')
    end
catch
    warndlg('Missing ANTs library. All functions in Analysis Tab cannot be used.')
end
end