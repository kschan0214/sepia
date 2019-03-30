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
%
%
function sepia
clear 

sepia_addpath;

global h

% set GUI window size
screenSize = get(0,'ScreenSize');
posLeft = round(screenSize(3)/4);
posBottom = round(screenSize(4)/6);
guiSizeHori = round(screenSize(3)/3);
guiSizeVert = round(screenSize(4)*2/3);
if guiSizeHori < 1000
    guiSizeHori = 1000;
end
if guiSizeVert < 700
    guiSizeVert = 700;
end

% create GUI 
h.fig=figure('Units','pixels','position',[posLeft posBottom guiSizeHori guiSizeVert],...
    'MenuBar','None','Toolbar','None','Name','Sepia GUI (v0.7.0)','NumberTitle','off');

% create Tabs for GUI
h.TabGroup          = uitabgroup(h.fig,'position',[.01 .01 0.98 0.98]);
h.Tabs.Sepia       = uitab(h.TabGroup,'Title','Sepia');
h.Tabs.phaseUnwrap  = uitab(h.TabGroup,'Title','Phase unwrapping');
h.Tabs.bkgRemoval   = uitab(h.TabGroup,'Title','Background field removal');
h.Tabs.qsm          = uitab(h.TabGroup,'Title','QSM');
h.Tabs.utility      = uitab(h.TabGroup,'Title','Utility');

% construct all tabs
%% Phase unwrapping tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.phaseUnwrap,              h,[0.01 0.8]);
% phase unwrap
h = qsmhub_handle_panel_phaseUnwrap(h.Tabs.phaseUnwrap,         h,[0.01 0.59]);

%% background field removal tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.bkgRemoval,               h,[0.01 0.8]);
% background field
h = qsmhub_handle_panel_bkgRemoval(h.Tabs.bkgRemoval,           h,[0.01 0.54]);

%% qsm tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.qsm,                      h,[0.01 0.8]);
% QSM
h = qsmhub_handle_panel_qsm(h.Tabs.qsm,                         h,[0.01 0.54]);

%% utility tab
h = qsmhub_handle_panel_Utility(h.Tabs.utility,                 h,[0.01 0.39]);

%% GUI with QSM one-stop station tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.Sepia,                   h,[0.01 0.8]);
% phase unwrap
h = qsmhub_handle_panel_phaseUnwrap(h.Tabs.Sepia,              h,[0.01 0.59]);
% background field
h = qsmhub_handle_panel_bkgRemoval(h.Tabs.Sepia,               h,[0.01 0.33]);
% QSM
h = qsmhub_handle_panel_qsm(h.Tabs.Sepia,                      h,[0.01 0.07]);

% Start button
h.pushbutton_start = uicontrol('Parent',h.Tabs.Sepia,...
    'Style','pushbutton',...
    'String','Start',...
    'units','normalized','Position',[0.85 0.01 0.1 0.05],...
    'backgroundcolor',get(h.fig,'color'));
% GPU checkbox
h.checkbox_gpu = uicontrol('Parent',h.Tabs.Sepia,...
    'Style','checkbox',...
    'String','Enable GPU computation',...
    'units','normalized','Position',[0.01 0.01 0.4 0.05],...
    'backgroundcolor',get(h.fig,'color'), 'Enable','off',...
    'TooltipString',['Enable to use GPU for some of the algorithms in sepia. ' ...
                     'Your GPU has to be detectable in Matlab in order to use this feature.']);
if gpuDeviceCount > 0
    set(h.checkbox_gpu, 'Enable', 'on');
end

%% Set Callback functions
set(h.TabGroup,      	'SelectionChangedFcn', {@SwitchTab_Callback})
set(h.pushbutton_start,	'Callback',            {@PushbuttonStart_Callback});

end

%% Callback functions
function SwitchTab_Callback(source,eventdata)
% switch parent handle of StepsPanel based on current tab

global h

switch eventdata.NewValue.Title
    
    % QSM one-stop station tab
    case 'Sepia'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.Sepia);
            % This tab supports both DICOM and NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains DICOM (both magnitude and phase files under the same directory) or NIfTI (*phase*.nii* and *magn*.nii*) files');
            % BET is supported with this tab
            set(h.dataIO.checkbox.brainExtraction,'Enable','on');
            % phase invert is supported with this tab
            set(h.dataIO.checkbox.invertPhase,'Enable','on');
            % input data 1
            set(h.dataIO.text.inputData1,'String','or Phase data:');
%             set(h.dataIO.edit.inputData1,'Enable','off');
            % input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
            set(h.dataIO.edit.inputData2,'Enable','on');
            set(h.dataIO.button.inputData2,'Enable','on');
            % input data 3
            set(h.dataIO.text.inputData3,'String','Weights:');
            set(h.dataIO.edit.inputData3,'Enable','on');
%             set(h.dataIO.edit.inputData3,'String',[]);
            set(h.dataIO.button.inputData3,'Enable','on');
            
        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.Sepia,'Position',[0.01 0.59 0.95 0.2]);
        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.Sepia,'Position',[0.01 0.33 0.95 0.25]);
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.Sepia,'Position',[0.01 0.07 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.Sepia);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.Sepia);
%         set(h.checkbox_gpu,             'Enable', 'on');
        
    % Phase unwrapping tab
    case 'Phase unwrapping'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.phaseUnwrap);
            % This tab supports both DICOM and NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains DICOM (both magnitude and phase files under the same directory) or NIfTI (*phase*.nii* and *magn*.nii*) files');
            % BET is supported with this tab
            set(h.dataIO.checkbox.brainExtraction,'Enable','on');
            % phase invert is supported with this tab
            set(h.dataIO.checkbox.invertPhase,'Enable','on');
            % input data 1
            set(h.dataIO.text.inputData1,'String','or Phase data:');
            % input data 2
%           set(h.dataIO.text.inputData2,'String','Magn. data:');
            set(h.dataIO.edit.inputData2,'Enable','on');
            set(h.dataIO.button.inputData2,'Enable','on');
            % input data 3
            set(h.dataIO.text.inputData3,'String','Weights:');
            set(h.dataIO.edit.inputData3,'Enable','off');
            set(h.dataIO.button.inputData3,'Enable','off');
            set(h.dataIO.edit.inputData3,'String',[]);

        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.phaseUnwrap,'Position',[0.01 0.59 0.95 0.2]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.phaseUnwrap);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.phaseUnwrap);
%         set(h.checkbox_gpu,             'Enable', 'off');
%         set(h.checkbox_gpu,             'Value',0);
        
    % background field removal tab    
    case 'Background field removal'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.bkgRemoval);
            % This tab supports only NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains the unwrapped total field map NIfTI (*total-field*.nii*) files');
            % no BET support with this tab
            set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            set(h.dataIO.edit.maskdir,              'Enable','on');
            set(h.dataIO.button.maskdir,            'Enable','on');
            % phase invert is not supported with this tab
            set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
            % input data 1
            set(h.dataIO.text.inputData1,'String','or Total field:');
            % input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
            set(h.dataIO.edit.inputData2,'Enable','off');
            set(h.dataIO.button.inputData2,'Enable','off');
            set(h.dataIO.edit.inputData2,'String',[]);
            % input data 3
            set(h.dataIO.text.inputData3,'String','Noise SD:');
            set(h.dataIO.edit.inputData3,'Enable','on');
            set(h.dataIO.button.inputData3,'Enable','on');

        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.bkgRemoval,'Position',[0.01 0.54 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.bkgRemoval);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.bkgRemoval);
%         set(h.checkbox_gpu,             'Enable', 'on');

    % qsm tab    
    case 'QSM'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.qsm);
            % This tab supports only NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains the local field map NIfTI (*local-field*.nii*) files; for some QSM methods additional file(s) may also be needed (e.g. *magn*.nii* and *weights*.nii*)');
            % no BET support with this tab
            set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            set(h.dataIO.edit.maskdir,              'Enable','on');
            set(h.dataIO.button.maskdir,            'Enable','on');
            % phase invert is not supported with this tab
            set(h.dataIO.checkbox.invertPhase,      'Enable','off','Value',0);
            % input data 1
            set(h.dataIO.text.inputData1,'String','or Local field:');
            % input data 2
%             set(h.dataIO.text.inputData2,'String','Magn. data:');
            set(h.dataIO.edit.inputData2,'Enable','on');
            set(h.dataIO.button.inputData2,'Enable','on');
            % input data 3
            set(h.dataIO.text.inputData3,'String','Weights:');
            set(h.dataIO.edit.inputData3,'Enable','on');
            set(h.dataIO.button.inputData3,'Enable','on');
%             % input header
%             set(h.dataIO.text.inputHeader,'String','Header:');
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.qsm,'Position',[0.01 0.54 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.qsm);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.qsm);
%         set(h.checkbox_gpu,             'Enable', 'on');
        
end

end

function PushbuttonStart_Callback(source,eventdata)
% core of Sepia GUI is to create a .m script to execute a command-based function

global h

% Disable the pushbutton to prevent double clicks
set(source,'Enable','off');

% get the tab name of the standalone
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
output_index = strfind(outputBasename, filesep);
outputDir = outputBasename(1:output_index(end));
% if the output directory does not exist then create the directory
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

% create a new m file
fid = fopen([outputDir filesep 'sepia_log.m'],'w');

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
fprintf(fid,'algorParam.general.isBET = %i ;\n'     ,get(h.dataIO.checkbox.brainExtraction, 'Value'));
fprintf(fid,'algorParam.general.isInvert = %i ;\n'  ,get(h.dataIO.checkbox.invertPhase,     'Value'));
fprintf(fid,'algorParam.general.isGPU = %i ;\n'     ,get(h.checkbox_gpu,                    'Value'));

% phase unwrap algorithm parameters
if strcmpi(tab,'Sepia') || strcmpi(tab,'Phase unwrapping')
    fprintf(fid,'%% Phase unwrapping algorithm parameters\n');
    % echo phase combine method
    fprintf(fid,'algorParam.unwrap.echoCombMethod = ''%s'' ;\n'     ,h.phaseUnwrap.popup.phaseCombMethod.String{h.phaseUnwrap.popup.phaseCombMethod.Value,1});
    % unwrap method
    switch h.phaseUnwrap.popup.phaseUnwrap.String{h.phaseUnwrap.popup.phaseUnwrap.Value,1}
        case 'Region growing'
            fprintf(fid,'algorParam.unwrap.unwrapMethod = ''%s'' ;\n'     ,'rg');

        case 'Graphcut'
            fprintf(fid,'algorParam.unwrap.unwrapMethod = ''%s'' ;\n'     ,'gc');
            fprintf(fid,'algorParam.unwrap.subsampling = %i ;\n'          ,1);

        case 'Laplacian STI suite'
            fprintf(fid,'algorParam.unwrap.unwrapMethod = ''%s'' ;\n'     ,'laplacian_stisuite');

        case '3D best path'
            fprintf(fid,'algorParam.unwrap.unwrapMethod = ''%s'' ;\n'     ,'bestpath3d');
    end
    % eddy current correction
    fprintf(fid,'algorParam.unwrap.isEddyCorrect = %i ;\n'     ,get(h.phaseUnwrap.checkbox.eddyCorrect,'Value'));
    % exclusion mask threshold
    if get(h.phaseUnwrap.checkbox.excludeMask,'Value')
        fprintf(fid,'algorParam.unwrap.excludeMaskThreshold = %g ;\n'     ,str2double(get(h.phaseUnwrap.edit.excludeMask,'String')));
    else
        fprintf(fid,'algorParam.unwrap.excludeMaskThreshold = Inf ;\n');
    end
end
    
% background field removal algorithm parameters
if strcmpi(tab,'Sepia') || strcmpi(tab,'Background field removal')
    fprintf(fid,'%% Background field removal algorithm parameters\n');
    % polyfit
    fprintf(fid,'algorParam.bfr.refine = %i ;\n'    ,get(h.bkgRemoval.checkbox.bkgRemoval,'Value'));
    % set parameters for selected method
    switch h.bkgRemoval.popup.bkgRemoval.String{h.bkgRemoval.popup.bkgRemoval.Value,1}
        case 'LBV'
            fprintf(fid,'algorParam.bfr.method = ''%s'' ;\n'     ,'lbv');
            fprintf(fid,'algorParam.bfr.tol = %s ;\n'	,(get(h.bkgRemoval.LBV.edit.tol,	'String')));
            fprintf(fid,'algorParam.bfr.depth = %s ;\n'	,(get(h.bkgRemoval.LBV.edit.depth,  'String')));
            fprintf(fid,'algorParam.bfr.peel = %s ;\n' 	,(get(h.bkgRemoval.LBV.edit.peel,	'String')));

        case 'PDF'
            fprintf(fid,'algorParam.bfr.method = ''%s'' ;\n'     ,'pdf');
            fprintf(fid,'algorParam.bfr.tol = %s ;\n'       ,(get(h.bkgRemoval.PDF.edit.tol,  	'String')));
            fprintf(fid,'algorParam.bfr.iteration = %s ;\n'	,(get(h.bkgRemoval.PDF.edit.maxIter,'String')));
            fprintf(fid,'algorParam.bfr.padSize = %s ;\n' 	,(get(h.bkgRemoval.PDF.edit.padSize,'String')));

        case 'RESHARP'
            fprintf(fid,'algorParam.bfr.method = ''%s'' ;\n'     ,'resharp');
            fprintf(fid,'algorParam.bfr.radius = %s ;\n'	,(get(h.bkgRemoval.RESHARP.edit.radius,  'String')));
            fprintf(fid,'algorParam.bfr.alpha = %s ;\n' 	,(get(h.bkgRemoval.RESHARP.edit.lambda,  'String')));

        case 'SHARP'
            fprintf(fid,'algorParam.bfr.method = ''%s'' ;\n'    ,'sharp');
            fprintf(fid,'algorParam.bfr.radius = %s ;\n'	,(get(h.bkgRemoval.SHARP.edit.radius,    'String')));
            fprintf(fid,'algorParam.bfr.threshold = %s ;\n'	,(get(h.bkgRemoval.SHARP.edit.threshold, 'String')));

        case 'VSHARP'
            fprintf(fid,'algorParam.bfr.method = ''%s'' ;\n'    ,'vsharp');
            fprintf(fid,'algorParam.bfr.radius = [%s:-1:%s] ;\n' ,(get(h.bkgRemoval.VSHARP.edit.maxRadius,'String')),...
                                                                   get(h.bkgRemoval.VSHARP.edit.minRadius,'String'));

        case 'iHARPERELLA'
            fprintf(fid,'algorParam.bfr.method = ''%s'' ;\n'    ,'iharperella');
            fprintf(fid,'algorParam.bfr.iteration = %s ;\n' ,(get(h.bkgRemoval.iHARPERELLA.edit.maxIter,'String')));

        case 'VSHARP STI suite'
            fprintf(fid,'algorParam.bfr.method = ''%s'' ;\n'    ,'vsharpsti');
            fprintf(fid,'algorParam.bfr.radius = %s ;\n' ,(get(h.bkgRemoval.VSHARPSTI.edit.smvSize,'String')));

    end
end

% QSM algorithm parameters
if strcmpi(tab,'Sepia') || strcmpi(tab,'QSM')
    fprintf(fid,'%% QSM algorithm parameters\n');
    % set parameters for selected method
    switch h.qsm.popup.qsm.String{h.qsm.popup.qsm.Value,1}
        case 'TKD'
            fprintf(fid,'algorParam.qsm.method = ''%s'' ;\n'	,'tkd');
            fprintf(fid,'algorParam.qsm.threshold = %s ;\n'     ,get(h.qsm.TKD.edit.threshold,'String'));
        
        case 'Closed-form solution'
            fprintf(fid,'algorParam.qsm.method = ''%s'' ;\n'    ,'closedforml2');
            fprintf(fid,'algorParam.qsm.lambda = %s ;\n'        ,get(h.qsm.cfs.edit.lambda,     'String'));
            fprintf(fid,'algorParam.qsm.optimise = %i ;\n'      ,get(h.qsm.cfs.checkbox.lambda, 'Value'));

        case 'STI suite iLSQR'
            fprintf(fid,'algorParam.qsm.method = ''%s'' ;\n'    ,'stisuiteilsqr');
            fprintf(fid,'algorParam.qsm.threshold = %s ;\n'         ,get(h.qsm.STIiLSQR.edit.threshold, 'String'));
            fprintf(fid,'algorParam.qsm.maxiter = %s ;\n'           ,get(h.qsm.STIiLSQR.edit.maxIter,   'String'));
            fprintf(fid,'algorParam.qsm.tol1 = %s ;\n'              ,get(h.qsm.STIiLSQR.edit.tol1,   'String'));
            fprintf(fid,'algorParam.qsm.tol2 = %s ;\n'              ,get(h.qsm.STIiLSQR.edit.tol2,      'String'));
            fprintf(fid,'algorParam.qsm.padsize = ones(1,3)*%s ;\n'	,get(h.qsm.STIiLSQR.edit.padSize,   'String'));

        case 'iLSQR'
            fprintf(fid,'algorParam.qsm.method = ''%s'' ;\n'    ,'ilsqr');
            fprintf(fid,'algorParam.qsm.tol = %s ;\n'           ,get(h.qsm.iLSQR.edit.tol,          'String'));
            fprintf(fid,'algorParam.qsm.maxiter = %s ;\n'      	,get(h.qsm.iLSQR.edit.maxIter,      'String'));
            fprintf(fid,'algorParam.qsm.lambda = %s ;\n'      	,get(h.qsm.iLSQR.edit.lambda,       'String'));
            fprintf(fid,'algorParam.qsm.optimise = %i ;\n'  	,get(h.qsm.iLSQR.checkbox.lambda,   'Value'));

        case 'FANSI'
            fprintf(fid,'algorParam.qsm.method = ''%s'' ;\n'    ,'fansi');
            fprintf(fid,'algorParam.qsm.tol = %s ;\n'           ,get(h.qsm.FANSI.edit.tol,      'String'));
            fprintf(fid,'algorParam.qsm.maxiter = %s ;\n'      	,get(h.qsm.FANSI.edit.maxIter,  'String'));
            fprintf(fid,'algorParam.qsm.lambda = %s ;\n'      	,get(h.qsm.FANSI.edit.lambda,   'String'));
            fprintf(fid,'algorParam.qsm.mu1 = %s ;\n'           ,get(h.qsm.FANSI.edit.mu,       'String'));
            fprintf(fid,'algorParam.qsm.mu2 = %s ;\n'           ,get(h.qsm.FANSI.edit.mu2,      'String'));
            fprintf(fid,'algorParam.qsm.solver = ''%s'' ;\n'  	,h.qsm.FANSI.popup.solver.String{h.qsm.FANSI.popup.solver.Value,1});
            fprintf(fid,'algorParam.qsm.constraint = ''%s'' ;\n',h.qsm.FANSI.popup.constraints.String{h.qsm.FANSI.popup.constraints.Value,1});

        case 'Star-QSM'
            fprintf(fid,'algorParam.qsm.method = ''%s'' ;\n'    ,'star');
            fprintf(fid,'algorParam.qsm.padsize = ones(1,3)*%s ;\n'	,get(h.qsm.Star.edit.padSize,   'String'));

        case 'MEDI'
            fprintf(fid,'algorParam.qsm.method = ''%s'' ;\n'    ,'medi_l1');
            fprintf(fid,'algorParam.qsm.lambda = %s ;\n'      	,get(h.qsm.MEDI.edit.lambda,        'String'));
            fprintf(fid,'algorParam.qsm.wData = %s ;\n'        	,get(h.qsm.MEDI.edit.weightData,    'String'));
            fprintf(fid,'algorParam.qsm.wGradient = %s ;\n'    	,get(h.qsm.MEDI.edit.weightGradient,'String'));
            fprintf(fid,'algorParam.qsm.zeropad = %s ;\n'      	,get(h.qsm.MEDI.edit.zeropad,       'String'));
            fprintf(fid,'algorParam.qsm.radius = %s ;\n'      	,get(h.qsm.MEDI.edit.smv_radius,    'String'));
            fprintf(fid,'algorParam.qsm.isSMV = %i ;\n'         ,get(h.qsm.MEDI.checkbox.smv,       'Value'));
            fprintf(fid,'algorParam.qsm.merit = %i ;\n'         ,get(h.qsm.MEDI.checkbox.merit,     'Value'));
            fprintf(fid,'algorParam.qsm.isLambdaCSF = %i ;\n'  	,get(h.qsm.MEDI.checkbox.lambda_csf,'Value'));
            fprintf(fid,'algorParam.qsm.lambdaCSF = %s ;\n'     ,get(h.qsm.MEDI.edit.lambda_csf,    'String'));  

    end
end

try
    switch tab
        case 'Sepia'
            fprintf(fid,'\nSepiaIOWrapper(input,output_basename,mask_filename,algorParam);\n');
            
        case 'Phase unwrapping'
            fprintf(fid,'\nUnwrapPhaseMacroIOWrapper(input,output_basename,mask_filename,algorParam);\n');
        
        case 'Background field removal'
            fprintf(fid,'\nBackgroundRemovalMacroIOWrapper(input,output_basename,mask_filename,algorParam);\n');
        
        case 'QSM'
            fprintf(fid,'\nQSMMacroIOWrapper(input,output_basename,mask_filename,algorParam);\n');
            
    end
    
    fclose(fid);
    
    % run process
    run([outputDir filesep 'sepia_log.m']);
    
    % re-enable the pushbutton
    set(source,'Enable','on');

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');

    rethrow(ME);
end


end
