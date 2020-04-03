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
% overloading. Nebertheless the GUI will automatically do this for you.
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
% Date last modified: 1 June 2018
% Date last modified: 2 March 2020 (v0.8.0)
%
function sepia 

% clear previous handles
clear global h

% make sure nothing is logged at the moment
diary off

% add path
sepia_addpath;

global h

%% create basic GUI
% set GUI window size
screenSize  = get(0,'ScreenSize');
posLeft     = round(screenSize(3)/4);
posBottom   = round(screenSize(4)/6);
guiSizeHori = round(screenSize(3)/3);
guiSizeVert = round(screenSize(4)*2/3);
if guiSizeHori < 1000
    guiSizeHori = 1000;
end
if guiSizeVert < 700
    guiSizeVert = 700;
end

% create GUI figure
h.fig=figure('Units','pixels','position',[posLeft posBottom guiSizeHori guiSizeVert],...
    'MenuBar','None','Toolbar','None','Name','SEPIA GUI (v0.8.0)','NumberTitle','off');

%% construct panels for each tab 
% create Tabs for GUI
h.TabGroup          = uitabgroup(h.fig,'position',[.01 .01 0.98 0.98]);
h.Tabs.Sepia        = uitab(h.TabGroup,'Title','SEPIA');
h.Tabs.phaseUnwrap  = uitab(h.TabGroup,'Title','Phase unwrapping');
h.Tabs.bkgRemoval   = uitab(h.TabGroup,'Title','Background field removal');
h.Tabs.qsm          = uitab(h.TabGroup,'Title','QSM');
h.Tabs.swismwi      = uitab(h.TabGroup,'Title','SWI/SMWI');
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
h = sepia_handle_panel_swi_dataIO(parent_curr,	h,[0.01 0.8]);
% Method
h = sepia_handle_panel_swi(parent_curr,       	h,[0.01 0.44]);

%% utility tab
h = sepia_handle_panel_Utility(h.Tabs.utility,	h,[0.01 0.39]);

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

% GPU checkbox
h.checkbox_gpu = uicontrol('Parent',h.Tabs.Sepia,...
    'Style','checkbox',...
    'String','Enable GPU computation',...
    'units','normalized','Position',[0.01 0.01 0.4 0.05],...
    'backgroundcolor',get(h.fig,'color'), ...
    'Enable','off','Visible','off',...
    'TooltipString',['Enable to use GPU for some of the algorithms in SEPIA. ' ...
                     'Your GPU has to be detectable in Matlab in order to use this feature.']);
if gpuDeviceCount > 0
    set(h.checkbox_gpu, 'Enable', 'on');
end

%% Set Callback functions
set(h.TabGroup,                 'SelectionChangedFcn', {@SwitchTab_Callback})
set(h.pushbutton_start,         'Callback',            {@PushbuttonStart_Callback});
set(h.pushbutton_loadConfig,    'Callback',            {@PushbuttonLoadConfig_Callback});

end

%% Callback functions

%% switching tabs
function SwitchTab_Callback(source,eventdata)
% switch parent handle of StepsPanel based on current tab

global h

% global uicontrol for all tabs
universial_handle = {h.StepsPanel.dataIO,...
                     h.pushbutton_start,...         % Start pushbutton
                     h.pushbutton_loadConfig,...    % load config pushbutton
                     h.checkbox_gpu};               % GPU checkbox

% Tab specific strings and tooltips
tooltip.input_dir{1} = 'Directory contains phase (*ph*.nii*), magnitude (*mag*.nii) & header (*header*.mat) (& mask, *mask*nii*) files';
tooltip.input_dir{2} = 'Directory contains the total field map (*total-field*.nii*) and SEPIA header (*header*.mat)';
tooltip.input_dir{3} = 'Directory contains the local field map (*local-field*.nii*) (depending on QSM algorithm, additional file(s) may also be needed, e.g. *mag*.nii* and *weights*.nii*)';

fieldString.inputData1{1}= 'or Phase:';
fieldString.inputData1{2}= 'or Total field:';
fieldString.inputData1{3}= 'or Local field:';

fieldString.inputData3{1}= '    Weights:';
fieldString.inputData3{2}= '    Noise SD:';

% change universal elements' parent except 'SWI/SMWI' tab and 'Utility' tab
if ~strcmpi(eventdata.NewValue.Title,'SWI/SMWI') && ~strcmpi(eventdata.NewValue.Title,'Utility')
    for k = 1:length(universial_handle)
        set(universial_handle{k}, 'Parent', source.SelectedTab);
    end
end

% Specify tab-specific content
switch eventdata.NewValue.Title
    
    % QSM one-stop station tab
    case 'SEPIA'
        % I/O
        % Change essential files if input is a directory
        set(h.dataIO.text.input,                'Tooltip',tooltip.input_dir{1});
        % BET is supported with this tab
        set(h.dataIO.checkbox.brainExtraction,  'Enable','on');
            % trigger followup callback to switch method panel
            feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
        % phase invert is supported with this tab
        set(h.dataIO.checkbox.invertPhase,      'Enable','on');
        % input data 1
        set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{1});
        % input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
        set(h.dataIO.edit.inputData2,           'Enable','on');
        set(h.dataIO.button.inputData2,         'Enable','on');
        % input data 3
        set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{1});
        set(h.dataIO.edit.inputData3,           'Enable','on');
        set(h.dataIO.button.inputData3,         'Enable','on');
            
        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.Sepia,'Position',[0.01 0.59 0.95 0.2]);
        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.Sepia,'Position',[0.01 0.33 0.95 0.25]);
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.Sepia,'Position',[0.01 0.07 0.95 0.25]);
        
    % Phase unwrapping tab
    case 'Phase unwrapping'
        % I/O
        % This tab supports both DICOM and NIfTI files
        set(h.dataIO.text.input,                'Tooltip',tooltip.input_dir{1});
        % BET is supported with this tab
        set(h.dataIO.checkbox.brainExtraction,  'Enable','on');
            % trigger followup callback to switch method panel
            feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
        % phase invert is supported with this tab
        set(h.dataIO.checkbox.invertPhase,      'Enable','on');
        % input data 1
        set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{1});
        % input data 2
%           set(h.dataIO.text.inputData2,'String','Magn. data:');
        set(h.dataIO.edit.inputData2,           'Enable','on');
        set(h.dataIO.button.inputData2,         'Enable','on');
        % input data 3
        set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{1});
        set(h.dataIO.edit.inputData3,           'Enable','off','String',[]);
        set(h.dataIO.button.inputData3,         'Enable','off');

        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.phaseUnwrap,'Position',[0.01 0.59 0.95 0.2]);
        
    % background field removal tab    
    case 'Background field removal'
        % I/O
        % This tab supports only NIfTI files
        set(h.dataIO.text.input,                'Tooltip',tooltip.input_dir{2});
        % no BET support with this tab
        set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            % trigger followup callback to switch method panel
            feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
        set(h.dataIO.edit.maskdir,              'Enable','on');
        set(h.dataIO.button.maskdir,            'Enable','on');
        % phase invert is not supported with this tab
        set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
        % input data 1
        set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{2});
        % input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
        set(h.dataIO.edit.inputData2,           'Enable','off','String',[]);
        set(h.dataIO.button.inputData2,         'Enable','off');
        % input data 3
        set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{2});
        set(h.dataIO.edit.inputData3,           'Enable','on');
        set(h.dataIO.button.inputData3,         'Enable','on');

        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.bkgRemoval,'Position',[0.01 0.54 0.95 0.25]);


    % qsm tab    
    case 'QSM'
        % I/O
        % This tab supports only NIfTI files
        set(h.dataIO.text.input,                'Tooltip',tooltip.input_dir{3});
        % no BET support with this tab
        set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            % trigger followup callback to switch method panel
            feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
        set(h.dataIO.edit.maskdir,              'Enable','on');
        set(h.dataIO.button.maskdir,            'Enable','on');
        % phase invert is not supported with this tab
        set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
        % input data 1
        set(h.dataIO.text.inputData1,           'String',fieldString.inputData1{3});
        % input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
        set(h.dataIO.edit.inputData2,           'Enable','on');
        set(h.dataIO.button.inputData2,         'Enable','on');
        % input data 3
        set(h.dataIO.text.inputData3,           'String',fieldString.inputData3{1});
        set(h.dataIO.edit.inputData3,           'Enable','on');
        set(h.dataIO.button.inputData3,         'Enable','on');
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.qsm,'Position',[0.01 0.54 0.95 0.25]);
        
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

% create a new m file
configFilename = [outputDir filesep 'sepia_config.m'];
if exist(configFilename,'file') == 2
    counter = 1;
    while exist(configFilename,'file') == 2
        suffix = ['_' num2str(counter)];
        configFilename = [outputDir filesep 'sepia_config' suffix '.m'];
        counter = counter + 1;
    end
end
fid = fopen(configFilename,'w');

%%%%%%%%%%%% Step 2: Write config file %%%%%%%%%%%% 
% general path
fprintf(fid,'%% add general Path\n');
fprintf(fid,'sepia_addpath\n\n');

fprintf(fid,'%% Input/Output filenames\n');
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
fprintf(fid,'algorParam.general.isBET       = %i ;\n'	,get(h.dataIO.checkbox.brainExtraction, 'Value'));
if get(h.dataIO.checkbox.brainExtraction, 'Value')
    fprintf(fid,'algorParam.general.fractional_threshold = %s ;\n'	,get(h.dataIO.edit.fractionalThres, 'String'));
    fprintf(fid,'algorParam.general.gradient_threshold   = %s ;\n'	,get(h.dataIO.edit.gradientThres,  	'String'));
end
fprintf(fid,'algorParam.general.isInvert    = %i ;\n'   ,get(h.dataIO.checkbox.invertPhase,     'Value'));
% fprintf(fid,'algorParam.general.isGPU       = %i ;\n'  	,get(h.checkbox_gpu,                    'Value'));

% phase unwrap algorithm parameters
if strcmpi(tab,'SEPIA') || strcmpi(tab,'Phase unwrapping')
    fprintf(fid,'%% Total field recovery algorithm parameters\n');
    % echo phase combine method
    fprintf(fid,'algorParam.unwrap.echoCombMethod = ''%s'' ;\n'     ,h.phaseUnwrap.popup.phaseCombMethod.String{h.phaseUnwrap.popup.phaseCombMethod.Value,1});
    % unwrap method
    fprintf(fid,'algorParam.unwrap.unwrapMethod   = ''%s'' ;\n'     ,h.phaseUnwrap.popup.phaseUnwrap.String{h.phaseUnwrap.popup.phaseUnwrap.Value,1});
    % eddy current correction
    fprintf(fid,'algorParam.unwrap.isEddyCorrect  = %i ;\n'         ,get(h.phaseUnwrap.checkbox.eddyCorrect,'Value'));
    % exclusion mask threshold
    if get(h.phaseUnwrap.checkbox.excludeMask,'Value')
        fprintf(fid,'algorParam.unwrap.excludeMaskThreshold = %g ;\n'     ,str2double(get(h.phaseUnwrap.edit.excludeMask,'String')));
        fprintf(fid,'algorParam.unwrap.excludeMethod        = ''%s'' ;\n' ,h.phaseUnwrap.popup.excludeMethod.String{h.phaseUnwrap.popup.excludeMethod.Value,1});
    end
    % save unwrapped echo phase
    fprintf(fid,'algorParam.unwrap.isSaveUnwrappedEcho = %i ;\n'     ,get(h.phaseUnwrap.checkbox.saveEchoPhase,'Value'));
end
                
% background field removal algorithm parameters
if strcmpi(tab,'SEPIA') || strcmpi(tab,'Background field removal')
    fprintf(fid,'%% Background field removal algorithm parameters\n');
    % polyfit
    fprintf(fid,'algorParam.bfr.refine = %i ;\n'        ,get(h.bkgRemoval.checkbox.bkgRemoval,'Value'));
    % Erode local field
    fprintf(fid,'algorParam.bfr.erode_radius = %s ;\n'	,get(h.bkgRemoval.edit.imerode,'String'));
    
    % set parameters for selected method
    bfr_method = h.bkgRemoval.popup.bkgRemoval.String{h.bkgRemoval.popup.bkgRemoval.Value,1};
    fprintf(fid,'algorParam.bfr.method = ''%s'' ;\n'     ,bfr_method);
    for k = 1:length(methodBFRName)
        if strcmpi(bfr_method,methodBFRName{k})
            feval(config_BFR_function{k},h,'set',fid);
        end
    end

end

% QSM algorithm parameters
if strcmpi(tab,'SEPIA') || strcmpi(tab,'QSM')
    fprintf(fid,'%% QSM algorithm parameters\n');
    
    fprintf(fid,'algorParam.qsm.reference_tissue = ''%s'' ;\n'   ,h.qsm.popup.tissue.String{h.qsm.popup.tissue.Value,1});
    
    qsm_method = h.qsm.popup.qsm.String{h.qsm.popup.qsm.Value,1};
    fprintf(fid,'algorParam.qsm.method = ''%s'' ;\n'     ,qsm_method);
    for k = 1:length(methodQSMName)
        if strcmpi(qsm_method,methodQSMName{k})
            feval(config_QSM_function{k},h,'set',fid);
        end
    end

end

% Determine application based on Tab
switch tab
    case 'SEPIA'
        fprintf(fid,'\nSepiaIOWrapper(input,output_basename,mask_filename,algorParam);\n');

    case 'Phase unwrapping'
        fprintf(fid,'\nUnwrapPhaseMacroIOWrapper(input,output_basename,mask_filename,algorParam);\n');

    case 'Background field removal'
        fprintf(fid,'\nBackgroundRemovalMacroIOWrapper(input,output_basename,mask_filename,algorParam);\n');

    case 'QSM'
        fprintf(fid,'\nQSMMacroIOWrapper(input,output_basename,mask_filename,algorParam);\n');

end

fclose(fid);

% log command window display to a text file
logFilename = [outputDir filesep 'run_sepia.log'];
if exist(logFilename,'file') == 2
    counter = 1;
    while exist(logFilename,'file') == 2
        suffix = ['_' num2str(counter)];
        logFilename = [outputDir filesep 'run_sepia' suffix '.log'];
        counter = counter + 1;
    end
end
diary(logFilename)

try
    % run process
    run(configFilename);
    
    % turn off the log
    diary off
    
    % re-enable the pushbutton
    set(source,'Enable','on');

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    
    % close log file
    disp('There was an error! Please check the command window/error message file for more information.');
    diary off
    
    % open a new text file for error message
    errorMessageFilename = [outputDir filesep 'run_sepia.error'];
    if exist(errorMessageFilename,'file') == 2
        counter = 1;
        while exist(errorMessageFilename,'file') == 2
            suffix = ['_' num2str(counter)];
            errorMessageFilename = [outputDir filesep 'run_sepia' suffix '.error'];
            counter = counter + 1;
        end
    end
    fid = fopen(errorMessageFilename,'w');
    fprintf(fid,'The identifier was:\n%s\n\n',ME.identifier);
    fprintf(fid,'The message was:\n\n');
    msgString = getReport(ME,'extended','hyperlinks','off');
    fprintf(fid,'%s',msgString);
    fclose(fid);
    
    % rethrow the error message to command window
    rethrow(ME);
    
end


end

function PushbuttonLoadConfig_Callback(source,eventdata)

global h

% only read m file
[config_filename,pathDir] = uigetfile({'*.m'},'Select a SEPIA config file');

% if file is specified then read it
if exist(fullfile(pathDir,config_filename),'file')
    set_config_Callback(fullfile(pathDir,config_filename),h);
end

tab = h.StepsPanel.dataIO.Parent.Title;

switch tab
    
    % QSM one-stop station tab
    case 'SEPIA'
        % I/O
        % BET is supported with this tab
        set(h.dataIO.checkbox.brainExtraction,  'Enable','on');
            % trigger followup callback to switch method panel
            feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
        % phase invert is supported with this tab
        set(h.dataIO.checkbox.invertPhase,      'Enable','on');
        % input data 1
        set(h.dataIO.edit.inputData2,           'Enable','on');
        set(h.dataIO.button.inputData2,         'Enable','on');
        % input data 3
        set(h.dataIO.edit.inputData3,           'Enable','on');
        set(h.dataIO.button.inputData3,         'Enable','on');
            
        
    % Phase unwrapping tab
    case 'Phase unwrapping'
        % I/O
        % BET is supported with this tab
        set(h.dataIO.checkbox.brainExtraction,  'Enable','on');
            % trigger followup callback to switch method panel
            feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
        % phase invert is supported with this tab
        set(h.dataIO.checkbox.invertPhase,      'Enable','on');
        % input data 2
%           set(h.dataIO.text.inputData2,'String','Magn. data:');
        set(h.dataIO.edit.inputData2,           'Enable','on');
        set(h.dataIO.button.inputData2,         'Enable','on');
        % input data 3
        set(h.dataIO.edit.inputData3,           'Enable','off','String',[]);
        set(h.dataIO.button.inputData3,         'Enable','off');

        
    % background field removal tab    
    case 'Background field removal'
        % no BET support with this tab
        set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            % trigger followup callback to switch method panel
            feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
        set(h.dataIO.edit.maskdir,              'Enable','on');
        set(h.dataIO.button.maskdir,            'Enable','on');
        % phase invert is not supported with this tab
        set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
        % input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
        set(h.dataIO.edit.inputData2,           'Enable','off','String',[]);
        set(h.dataIO.button.inputData2,         'Enable','off');
        % input data 3
        set(h.dataIO.edit.inputData3,           'Enable','on');
        set(h.dataIO.button.inputData3,         'Enable','on');

    % qsm tab    
    case 'QSM'
        
        % no BET support with this tab
        set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            % trigger followup callback to switch method panel
            feval(h.dataIO.checkbox.brainExtraction.Callback{1},h.dataIO.checkbox.brainExtraction,[],h);
        set(h.dataIO.edit.maskdir,              'Enable','on');
        set(h.dataIO.button.maskdir,            'Enable','on');
        % phase invert is not supported with this tab
        set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
       
        set(h.dataIO.edit.inputData2,           'Enable','on');
        set(h.dataIO.button.inputData2,         'Enable','on');
        % input data 3
        set(h.dataIO.edit.inputData3,           'Enable','on');
        set(h.dataIO.button.inputData3,         'Enable','on');
        
end

end