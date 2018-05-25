%% qsm_hub
%
% Description: This is a GUI of qsm_hub, which is a pipeline control tool
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
% Date last modified: 20 May 2018
%
%
function qsm_hub
clear 

% qsm_hub_AddPath
qsm_hub_AddMethodPath;

global h

% set GUI window size
screenSize = get(0,'ScreenSize');
posLeft = round(screenSize(3)/4);
posBottom = round(screenSize(4)/6);
guiSizeHori = round(screenSize(3)/3);
guiSizeVert = round(screenSize(4)*2/3);
if guiSizeHori < 950
    guiSizeHori = 950;
end
if guiSizeVert < 700
    guiSizeVert = 700;
end

% create GUI 
fig=figure('Units','pixels','position',[posLeft posBottom guiSizeHori guiSizeVert],...
    'MenuBar','None','Toolbar','None','Name','QSM hub (beta)','NumberTitle','off');

% create Tabs for GUI
h.TabGroup          = uitabgroup(fig,'position',[.01 .01 0.99 0.99]);
h.Tabs.QSMHub       = uitab(h.TabGroup,'Title','One-stop QSM processing');
h.Tabs.phaseUnwrap  = uitab(h.TabGroup,'Title','Phase unwrapping');
h.Tabs.bkgRemoval   = uitab(h.TabGroup,'Title','Background field removal');
h.Tabs.qsm          = uitab(h.TabGroup,'Title','QSM');
h.Tabs.utility      = uitab(h.TabGroup,'Title','Utility');

% initialise all tabs
%% Phase unwrapping tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.phaseUnwrap,fig,h,[0.01 0.8]);
% phase unwrap
h = qsmhub_handle_panel_phaseUnwrap(h.Tabs.phaseUnwrap,fig,h,[0.01 0.59]);

%% background field removal tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.bkgRemoval,fig,h,[0.01 0.8]);
% background field
h = qsmhub_handle_panel_bkgRemoval(h.Tabs.bkgRemoval,fig,h,[0.01 0.54]);

%% qsm tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.qsm,fig,h,[0.01 0.8]);
% QSM
h = qsmhub_handle_panel_qsm(h.Tabs.qsm,fig,h,[0.01 0.54]);

%% utility tab
% get header
h = qsmhub_handle_panel_utility_get_header(h.Tabs.utility,fig,h,[0.01 0.7]);
% get ventricle mask
h = qsmhub_handle_panel_utility_mask_ventricle(h.Tabs.utility,fig,h,[0.01 0.45]);

%% GUI with QSM one-stop station tab
% I/O
h = qsmhub_handle_panel_dataIO(h.Tabs.QSMHub,fig,h,[0.01 0.8]);
% phase unwrap
h = qsmhub_handle_panel_phaseUnwrap(h.Tabs.QSMHub,fig,h,[0.01 0.59]);
% background field
h = qsmhub_handle_panel_bkgRemoval(h.Tabs.QSMHub,fig,h,[0.01 0.33]);
% QSM
h = qsmhub_handle_panel_qsm(h.Tabs.QSMHub,fig,h,[0.01 0.07]);
% Start button
h.pushbutton_start = uicontrol('Parent',h.Tabs.QSMHub,'Style','pushbutton',...
    'String','Start',...
    'units','normalized','Position',[0.85 0.01 0.1 0.05],...
    'backgroundcolor',get(fig,'color'));
% GPU checkbox
h.checkbox_gpu = uicontrol('Parent',h.Tabs.QSMHub,'Style','checkbox',...
    'String','Enable GPU computation',...
    'units','normalized','Position',[0.01 0.01 0.4 0.05],...
    'backgroundcolor',get(fig,'color'), 'Enable','off',...
    'TooltipString',['Enable to use GPU for some of the algorithms in qsm_hub. ' ...
                     'Your GPU has to be detectable in Matlab in order to use this feature.']);
if gpuDeviceCount > 0
    set(h.checkbox_gpu, 'Enable', 'on');
end

%% Set Callback functions
h = SetAllCallbacks(h);
end

%% utils functions
function h=SetAllCallbacks(h)
set(h.TabGroup,                             'SelectionChangedFcn',  {@SwitchTab_Callback})
set(h.dataIO.button.input,                  'Callback',             {@ButtonOpen_Callback,'input'});
set(h.dataIO.button.output,                 'Callback',             {@ButtonOpen_Callback,'output'});
set(h.dataIO.checkbox.brainExtraction,      'Callback',             {@CheckboxBrainExtraction_Callback});
set(h.dataIO.button.maskdir,                'Callback',             {@ButtonOpen_Callback,'mask'});
set(h.phaseUnwrap.checkbox.excludeMask,     'Callback',             {@CheckboxEditPair_Callback,h.phaseUnwrap.edit.excludeMask,1});
set(h.phaseUnwrap.edit.excludeMask,         'Callback',             {@EditInputMinMax_Callback,0,0,1});
set(h.bkgRemoval.popup.bkgRemoval,      	'Callback',             {@PopupBkgRemoval_Callback});
set(h.bkgRemoval.LBV.edit.depth,            'Callback',             {@EditInputMinMax_Callback,1,-1});
set(h.bkgRemoval.LBV.edit.peel,             'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.bkgRemoval.LBV.edit.tol,              'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.bkgRemoval.PDF.edit.maxIter,        	'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.bkgRemoval.PDF.edit.tol,              'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.bkgRemoval.PDF.edit.padSize,      	'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.bkgRemoval.RESHARP.edit.lambda,     	'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.bkgRemoval.RESHARP.edit.radius,     	'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.bkgRemoval.SHARP.edit.radius,     	'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.bkgRemoval.SHARP.edit.threshold,    	'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.bkgRemoval.VSHARP.edit.minRadius,   	'Callback',             {@EditVSHARPRadius_Callback});
set(h.bkgRemoval.VSHARP.edit.maxRadius,   	'Callback',             {@EditVSHARPRadius_Callback});
set(h.bkgRemoval.VSHARPSTI.edit.smvSize,    'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.bkgRemoval.iHARPERELLA.edit.maxIter,	'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.qsm.popup.qsm,                        'Callback',             {@PopupQSM_Callback});
set(h.qsm.cfs.checkbox.lambda,              'Callback',             {@CheckboxEditPair_Callback,h.qsm.cfs.edit.lambda,0});
set(h.qsm.iLSQR.checkbox.lambda,            'Callback',             {@CheckboxEditPair_Callback,h.qsm.iLSQR.edit.lambda,0});
set(h.qsm.TKD.edit.threshold,               'Callback',             {@EditInputMinMax_Callback,0,0,1});
set(h.qsm.cfs.edit.lambda,                  'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.Star.edit.padSize,                'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.qsm.iLSQR.edit.lambda,                'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.iLSQR.edit.maxIter,               'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.qsm.iLSQR.edit.tol,                   'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.STIiLSQR.edit.maxIter,            'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.qsm.STIiLSQR.edit.padSize,            'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.qsm.STIiLSQR.edit.threshold,          'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.STIiLSQR.edit.tol1,               'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.STIiLSQR.edit.tol2,               'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.FANSI.edit.lambda,                'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.FANSI.edit.maxIter,               'Callback',             {@EditInputMinMax_Callback,1,0});
set(h.qsm.FANSI.edit.mu,                    'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.FANSI.edit.mu2,                   'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.FANSI.edit.tol,                   'Callback',             {@EditInputMinMax_Callback,0,0});
set(h.qsm.MEDI.checkbox.smv,                'Callback',             {@CheckboxEditPair_Callback,h.qsm.MEDI.edit.smv_radius,1});
set(h.qsm.MEDI.checkbox.lambda_csf,         'Callback',             {@CheckboxEditPair_Callback,h.qsm.MEDI.edit.lambda_csf,1});
set(h.pushbutton_start,                     'Callback',             {@PushbuttonStart_Callback});

set(h.Utility.getHeader.button.teFile,      'Callback',             {@ButtonOpen_Utility_getHeader_Callback,'te'});
set(h.Utility.getHeader.button.dicomInput,  'Callback',             {@ButtonOpen_Utility_getHeader_Callback,'dicom'});
set(h.Utility.getHeader.button.niftiInput,  'Callback',             {@ButtonOpen_Utility_getHeader_Callback,'nitfi'});
set(h.Utility.getHeader.button.outputDir,   'Callback',             {@ButtonOpen_Utility_getHeader_Callback,'output'});
set(h.Utility.getHeader.button.run,         'Callback',             {@PushbuttonRun_Utility_getHeader_Callback});

set(h.Utility.csfMask.button.teFile,        'Callback',             {@ButtonOpen_Utility_csfMask_Callback,'te'});
set(h.Utility.csfMask.button.niftiInput,    'Callback',             {@ButtonOpen_Utility_csfMask_Callback,'nitfi'});
set(h.Utility.csfMask.button.maskInput,     'Callback',             {@ButtonOpen_Utility_csfMask_Callback,'mask'});
set(h.Utility.csfMask.button.outputDir,     'Callback',             {@ButtonOpen_Utility_csfMask_Callback,'output'});
set(h.Utility.csfMask.button.run,           'Callback',             {@PushbuttonRun_Utility_csfMask_Callback});
end

%% Callback functions
%% Common callback functions
function CheckboxEditPair_Callback(source,eventdata,handleToBeDisable,trueValue)
% enable/disable edit fields via checkbox
    % compare source value to trueValue
    % e.g. if an edit field needed to be disable by checking an checkbox
    % then trueValue of the checkbox is 0
    if source.Value == trueValue
        % if source is equal to trueValue then enables target handle
        set(handleToBeDisable,'Enable','on');
    else
        % if source do not equal to trueValue then disables target handle 
        set(handleToBeDisable,'Enable','off');
    end

end

function EditInputMinMax_Callback(source,eventdata,isIntegerInput,lb,ub)
% set the min/max input allowed in edit fields

    % check minimum
    if str2double(source.String)<lb
        source.String = num2str(lb);
    end
    
    % if maximum is setted then check maximum
    if nargin==5
        if str2double(source.String)>ub
            source.String = num2str(ub);
        end
    end
    
    % make sure the input is interger for some fields
    if isIntegerInput
        source.String = num2str(round(str2double(source.String)));
    end
    
end

%% Specific callback functions
function SwitchTab_Callback(source,eventdata)
% switch parent handle of StepsPanel based on current tab

global h

switch eventdata.NewValue.Title
    
    % QSM one-stop station tab
    case 'One-stop QSM processing'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.QSMHub);
            % This tab supports both DICOM and NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains DICOM (both magnitude and phase files under the same directory) or NIfTI (*phase*.nii* and *magn*.nii*) files');
            % BET is supported with this tab
            set(h.dataIO.checkbox.brainExtraction,'Enable','on');
        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.QSMHub,'Position',[0.01 0.59 0.95 0.2]);
        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.QSMHub,'Position',[0.01 0.33 0.95 0.25]);
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.QSMHub,'Position',[0.01 0.07 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.QSMHub);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.QSMHub);
        
    % Phase unwrapping tab
    case 'Phase unwrapping'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.phaseUnwrap);
            % This tab supports both DICOM and NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains DICOM (both magnitude and phase files under the same directory) or NIfTI (*phase*.nii* and *magn*.nii*) files');
            % BET is supported with this tab
            set(h.dataIO.checkbox.brainExtraction,'Enable','on');
        % phase unwrap
        set(h.StepsPanel.phaseUnwrap,   'Parent',h.Tabs.phaseUnwrap,'Position',[0.01 0.59 0.95 0.2]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.phaseUnwrap);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.phaseUnwrap);
        
    % background field removal tab    
    case 'Background field removal'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.bkgRemoval);
            % This tab supports only NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains the unwrapped total field map NIfTI (*totalField*.nii*) files');
            % no BET support with this tab
            set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            set(h.dataIO.edit.maskdir,              'Enable','on');
            set(h.dataIO.button.maskdir,            'Enable','on');
        % background field
        set(h.StepsPanel.bkgRemoval,    'Parent',h.Tabs.bkgRemoval,'Position',[0.01 0.54 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.bkgRemoval);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.bkgRemoval);

    % qsm tab    7
    case 'QSM'
        % I/O
        set(h.StepsPanel.dataIO,        'Parent',h.Tabs.qsm);
            % This tab supports only NIfTI files
            set(h.dataIO.text.input, 'Tooltip',...
                'Input directory contains the local field map NIfTI (*localField*.nii*) files; for some QSM methods additional file(s) may also be needed (e.g. *magn*.nii* and *fieldmapSD*.nii*)');
            % no BET support with this tab
            set(h.dataIO.checkbox.brainExtraction,  'Enable','off','Value',0);
            set(h.dataIO.edit.maskdir,              'Enable','on');
            set(h.dataIO.button.maskdir,            'Enable','on');
        % QSM
        set(h.StepsPanel.qsm,           'Parent',h.Tabs.qsm,'Position',[0.01 0.54 0.95 0.25]);
        % Start pushbutton
        set(h.pushbutton_start,         'Parent',h.Tabs.qsm);
        % GPU checkbox
        set(h.checkbox_gpu,             'Parent',h.Tabs.qsm);
        
end

end

function ButtonOpen_Callback(source,eventdata,field)
% get directory and display it on GUI

global h

switch field
    case 'mask'
        % only read NIfTI file for mask
        [maskfileName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select mask file');

        if pathDir ~= 0
            set(h.dataIO.edit.maskdir,'String',fullfile(pathDir,maskfileName));
        end
        
    case 'input'
        % get directory for NIfTI or DICOM files
        pathDir = uigetdir;

        if pathDir ~= 0
            % set input edit field for display
            set(h.dataIO.edit.input,    'String',pathDir);
            % automatically set default output field
            set(h.dataIO.edit.output,   'String',[pathDir filesep 'output']);
        end
        
    case 'output'
        
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.dataIO.edit.output,'String',pathDir);
        end
end

end

function CheckboxBrainExtraction_Callback(source,eventdata)
% if BET checkbox is checked then empty mask edit field and disable open
% pushbutton

global h

if ~h.dataIO.checkbox.brainExtraction.Value
    set(h.dataIO.button.maskdir,'Enable','on');
    set(h.dataIO.edit.maskdir,  'Enable','on');
else
    set(h.dataIO.button.maskdir,'Enable','off');
    set(h.dataIO.edit.maskdir,  'Enable','off');
    set(h.dataIO.edit.maskdir,  'String','');
end

end

function PopupBkgRemoval_Callback(source,eventdata)
% display corresponding background field removal method's panel

global h

% get selected background removal method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.bkgRemoval.panel);
for kf = 1:length(fields)
    set(h.bkgRemoval.panel.(fields{kf}),    'Visible','off');
end

% switch on target panel
switch method
    case 'LBV'
        set(h.bkgRemoval.panel.LBV,         'Visible','on');
        
    case 'PDF'
        set(h.bkgRemoval.panel.PDF,         'Visible','on');

    case 'RESHARP'
        set(h.bkgRemoval.panel.RESHARP,     'Visible','on');

    case 'SHARP'
        set(h.bkgRemoval.panel.SHARP,       'Visible','on');

    case 'VSHARP'
        set(h.bkgRemoval.panel.VSHARP,      'Visible','on');

    case 'VSHARP STI suite'
        set(h.bkgRemoval.panel.VSHARPSTI,   'Visible','on');

    case 'iHARPERELLA'
        set(h.bkgRemoval.panel.iHARPERELLA, 'Visible','on');
end

end

function PopupQSM_Callback(source,eventdata)
% display corresponding QSM method's panel

global h

% get selected QSM method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.qsm.panel);
for kf = 1:length(fields)
    set(h.qsm.panel.(fields{kf}),   'Visible','off');
end

% switch on target panel
switch method
    case 'TKD'
        set(h.qsm.panel.TKD,        'Visible','on');

    case 'Closed-form solution'
        set(h.qsm.panel.cfs,        'Visible','on');

    case 'STI suite iLSQR'
        set(h.qsm.panel.STIiLSQR,   'Visible','on');

    case 'iLSQR'
        set(h.qsm.panel.iLSQR,      'Visible','on');

    case 'FANSI'
        set(h.qsm.panel.FANSI,      'Visible','on');

    case 'Star'
        set(h.qsm.panel.Star,       'Visible','on');

    case 'MEDI'
        set(h.qsm.panel.MEDI,       'Visible','on');
end

end

function EditVSHARPRadius_Callback(source,eventdata)
% constraint the minimum of maximum radius is always larger then the
% minimum radius

global h

% check minimum of minimum radius input
if str2double(h.bkgRemoval.VSHARP.edit.minRadius.String)<0
    h.bkgRemoval.VSHARP.edit.minRadius.String = num2str(0);
end

% if the minimum radius is not integer then rounds it to interger
h.bkgRemoval.VSHARP.edit.minRadius.String = num2str(round(str2double(h.bkgRemoval.VSHARP.edit.minRadius.String)));

% ensure maximum radius is always larger then minimum radius
if str2double(h.bkgRemoval.VSHARP.edit.maxRadius.String) <= str2double(h.bkgRemoval.VSHARP.edit.minRadius.String)
    h.bkgRemoval.VSHARP.edit.maxRadius.String = num2str(str2double(h.bkgRemoval.VSHARP.edit.minRadius.String) +1);
end

% if the maximum radius is not integer then rounds it to interger
h.bkgRemoval.VSHARP.edit.maxRadius.String = num2str(round(str2double(h.bkgRemoval.VSHARP.edit.maxRadius.String)));

end

function PushbuttonStart_Callback(source,eventdata)
% determine which step(s) is going to process

global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

% initialise all possible parameters
subsampling=1;
BFR_tol=1e-4;BFR_depth=4;BFR_peel=2;BFR_iteration=50;
% BFR_CGdefault=true;
BFR_padSize = 40;        
        
BFR_radius=4;BFR_alpha=0.01;BFR_threshold=0.03;
QSM_threshold=0.15;QSM_lambda=0.13;QSM_optimise=false;
QSM_tol=1e-3;QSM_maxiter=50;QSM_tol1=0.01;QSM_tol2=0.001;QSM_padsize=[4,4,4];
QSM_mu1=5e-5;QSM_solver='linear';QSM_constraint='tv';QSM_mu2=1;
QSM_radius=5;QSM_zeropad=0;QSM_wData=1;QSM_wGradient=1;QSM_lambdaCSF=100;
QSM_isSMV=false;QSM_merit=false;QSM_isLambdaCSF=false;

% get GPU enable
isGPU = get(h.checkbox_gpu,'Value');

% get I/O GUI input
inputDir        = get(h.dataIO.edit.input,'String');
outputDir       = get(h.dataIO.edit.output,'String');
maskFullName    = get(h.dataIO.edit.maskdir,'String');
isBET           = get(h.dataIO.checkbox.brainExtraction,'Value');

% get phase unwrap GUI input
phaseUnwrap     = h.phaseUnwrap.popup.phaseUnwrap.String{h.phaseUnwrap.popup.phaseUnwrap.Value,1};
isEddyCorrect   = get(h.phaseUnwrap.checkbox.eddyCorrect,'Value');
if get(h.phaseUnwrap.checkbox.excludeMask,'Value')
    excludeMaskThreshold = str2double(get(h.phaseUnwrap.edit.excludeMask,'String'));
else
    excludeMaskThreshold = 1;
end

% get background field removal GUI input
BFR             = h.bkgRemoval.popup.bkgRemoval.String{h.bkgRemoval.popup.bkgRemoval.Value,1};
refine          = get(h.bkgRemoval.checkbox.bkgRemoval,'Value');

% get QSM GUI input
QSM_method      = h.qsm.popup.qsm.String{h.qsm.popup.qsm.Value,1};

% match the phase unwrapping GUI input to QSMHub input format
switch phaseUnwrap
    case 'Region growing'
        phaseUnwrap = 'rg';
        
    case 'Graphcut'
        phaseUnwrap = 'gc';
        
    case 'Laplacian STI suite'
        phaseUnwrap = 'laplacian_stisuite';
        
    case '3D best path'
        phaseUnwrap = 'bestpath3d';
        
end

% match the background field removal GUI input to QSMHub input format
% get specific backgroud field removal algorithm parameters
switch BFR
    case 'LBV'
        BFR='lbv';
        try BFR_tol         = str2double(get(h.bkgRemoval.LBV.edit.tol,'String'));              catch; BFR_tol=1e-4;        end
        try BFR_depth       = str2double(get(h.bkgRemoval.LBV.edit.depth,'String'));            catch; BFR_depth=4;         end
        try BFR_peel        = str2double(get(h.bkgRemoval.LBV.edit.peel,'String'));             catch; BFR_peel=4;          end
        
    case 'PDF'
        BFR='pdf';
        try BFR_tol         = str2double(get(h.bkgRemoval.PDF.edit.tol,'String'));              catch; BFR_tol=1e-2;        end
        try BFR_iteration   = str2double(get(h.bkgRemoval.PDF.edit.maxIter,'String'));          catch; BFR_iteration=50;    end
        try BFR_padSize     = str2double(get(h.bkgRemoval.PDF.edit.padSize,'String'));          catch; BFR_iteration=40;    end
%         try BFR_CGdefault = h.popup_PDF_cgSolver.String{h.popup_PDF_cgSolver.Value,1}; catch; BFR_CGdefault=true; end

    case 'RESHARP'
        BFR='resharp';
        try BFR_radius      = str2double(get(h.bkgRemoval.RESHARP.edit.radius,'String'));       catch; BFR_radius=4;        end
        try BFR_alpha       = str2double(get(h.bkgRemoval.RESHARP.edit.lambda,'String'));       catch; BFR_alpha=0.01;      end
        
    case 'SHARP'
        BFR='sharp';
        try BFR_radius      = str2double(get(h.bkgRemoval.SHARP.edit.radius,'String'));         catch; BFR_radius=4;        end
        try BFR_threshold   = str2double(get(h.bkgRemoval.SHARP.edit.threshold,'String'));      catch; BFR_threshold=0.03;	end
        
    case 'VSHARP'
        BFR='vsharp';
        try maxRadius       = str2double(get(h.bkgRemoval.VSHARP.edit.maxRadius,'String'));     catch; maxRadius=10;        end
        try minRadius       = str2double(get(h.bkgRemoval.VSHARP.edit.minRadius,'String'));     catch; minRadius=3;         end
        
        BFR_radius = maxRadius:-1:minRadius;
        
    case 'iHARPERELLA'
        BFR='iharperella';
        try BFR_iteration   = str2double(get(h.bkgRemoval.iHARPERELLA.edit.maxIter,'String'));  catch; BFR_iteration=100;   end 
        
    case 'VSHARP STI suite'
        BFR='vsharpsti';
        try BFR_radius      = str2double(get(h.bkgRemoval.VSHARPSTI.edit.smvSize,'String'));    catch; BFR_radius=12;       end
        
end

% get QSM algorithm parameters
switch QSM_method
    case 'TKD'
        QSM_method='tkd';
        try QSM_threshold   = str2double(get(h.qsm.TKD.edit.threshold,'String'));       catch; QSM_threshold=0.15;  end
        
    case 'Closed-form solution'
        QSM_method='closedforml2';
        try QSM_lambda      = str2double(get(h.qsm.cfs.edit.lambda,'String'));          catch; QSM_lambda=0.13;     end
        try QSM_optimise    = get(h.qsm.cfs.checkbox.lambda,'Value');                  	catch; QSM_optimise=false;  end
        
    case 'STI suite iLSQR'
        QSM_method='stisuiteilsqr';
        try QSM_threshold	= str2double(get(h.qsm.STIiLSQR.edit.threshold,'String'));  catch; QSM_threshold=0.01;  end
        try QSM_maxiter     = str2double(get(h.qsm.STIiLSQR.edit.maxIter,'String'));	catch; QSM_maxiter=100;     end
        try QSM_tol1        = str2double(get(h.qsm.STIiLSQR.edit.tol1,'String'));       catch; QSM_tol1=0.01;       end
        try QSM_tol2        = str2double(get(h.qsm.STIiLSQR.edit.tol2,'String'));       catch; QSM_tol2=0.001;      end
        try QSM_padsize     = str2double(get(h.qsm.STIiLSQR.edit.padSize,'String'));    catch; QSM_padsize=4;       end
        
        QSM_padsize = [QSM_padsize,QSM_padsize,QSM_padsize];
        
    case 'iLSQR'
        QSM_method='ilsqr';
        try QSM_tol         = str2double(get(h.qsm.iLSQR.edit.tol,'String'));           catch; QSM_tol=0.001;       end
        try QSM_maxiter     = str2double(get(h.qsm.iLSQR.edit.maxIter,'String'));       catch; QSM_maxiter=100;     end
        try QSM_lambda      = str2double(get(h.qsm.iLSQR.edit.lambda,'String'));        catch; QSM_lambda=0.13;     end
        try QSM_optimise    = get(h.qsm.iLSQR.checkbox.lambda,'Value');                     catch; QSM_optimise=false;  end 
        
    case 'FANSI'
        QSM_method='fansi';
        try QSM_tol         = str2double(get(h.qsm.FANSI.edit.tol,'String'));           catch; QSM_tol=1;           end
        try QSM_lambda      = str2double(get(h.qsm.FANSI.edit.lambda,'String'));        catch; QSM_lambda=3e-5;     end
        try QSM_mu1         = str2double(get(h.qsm.FANSI.edit.mu,'String'));            catch; QSM_mu1=5e-5;        end
        try QSM_mu2         = str2double(get(h.qsm.FANSI.edit.mu2,'String'));           catch; QSM_mu2=1;           end
        try QSM_maxiter     = str2double(get(h.qsm.FANSI.edit.maxIter,'String'));       catch; QSM_maxiter=50;      end
        try 
            QSM_solver      = h.qsm.FANSI.popup.solver.String{h.qsm.FANSI.popup.solver.Value,1}; 
        catch
            QSM_solver      ='linear'; 
        end 
        try 
            QSM_constraint  = h.qsm.FANSI.popup.constraints.String{h.qsm.FANSI.popup.constraints.Value,1}; 
        catch
            QSM_constraint  ='tv'; 
        end 
        
    case 'Star'
        QSM_method='star';
        try QSM_padsize   = str2double(get(h.qsm.Star.edit.padSize,'String'));        catch; QSM_padsize=4;       end
        QSM_padsize = [QSM_padsize,QSM_padsize,QSM_padsize];
        
    case 'MEDI'
        QSM_method='medi_l1';
        try QSM_lambda      = str2double(get(h.qsm.MEDI.edit.lambda,'String'));         catch; QSM_lambda=1000;     end 
        try QSM_wData       = str2double(get(h.qsm.MEDI.edit.weightData,'String'));     catch; QSM_wData=1;         end 
        try QSM_wGradient   = str2double(get(h.qsm.MEDI.edit.weightGradient,'String')); catch; QSM_wGradient=1;     end 
        try QSM_zeropad     = str2double(get(h.qsm.MEDI.edit.zeropad,'String'));        catch; QSM_zeropad=0;       end 
        try QSM_radius      = str2double(get(h.qsm.MEDI.edit.smv_radius,'String'));     catch; QSM_radius=5;        end 
        try QSM_isSMV       = get(h.qsm.MEDI.checkbox.smv,'Value');                     catch; QSM_isSMV=0;         end 
        try QSM_merit       = get(h.qsm.MEDI.checkbox.merit,'Value');                   catch; QSM_merit=0;         end 
        try QSM_isLambdaCSF = get(h.qsm.MEDI.checkbox.lambda_csf,'Value');              catch; QSM_isLambdaCSF=0;   end 
        try QSM_lambdaCSF   = str2double(get(h.qsm.MEDI.edit.lambda_csf,'String'));     catch; QSM_lambdaCSF=100;   end 
        
end

% run the selected processing step and delay the error (if any)
try 
    % get the parent of current panel, this determine which script is about to run
    switch h.StepsPanel.dataIO.Parent.Title
        case 'One-stop QSM processing'
            % core of QSM one-stop processing
            QSMHub(inputDir,outputDir,'FSLBet',isBET,'mask',maskFullName,'unwrap',phaseUnwrap,...
                'Subsampling',subsampling,'BFR',BFR,'refine',refine,'BFR_tol',BFR_tol,...
                'depth',BFR_depth,'peel',BFR_peel,'BFR_iteration',BFR_iteration,'BFR_padsize',BFR_padSize,...
                'BFR_radius',BFR_radius,'BFR_alpha',BFR_alpha,'BFR_threshold',BFR_threshold,...
                'QSM',QSM_method,'QSM_threshold',QSM_threshold,'QSM_lambda',QSM_lambda,'QSM_optimise',QSM_optimise,...
                'QSM_tol',QSM_tol,'QSM_iteration',QSM_maxiter,'QSM_tol1',QSM_tol1,'QSM_tol2',QSM_tol2,...
                'QSM_padsize',QSM_padsize,'QSM_mu',QSM_mu1,'QSM_mu2',QSM_mu2,QSM_solver,QSM_constraint,'exclude_threshold',excludeMaskThreshold,...
                'QSM_zeropad',QSM_zeropad,'QSM_wData',QSM_wData,'QSM_wGradient',QSM_wGradient,'QSM_radius',QSM_radius,...
                'QSM_isSMV',QSM_isSMV,'QSM_merit',QSM_merit,'QSM_isLambdaCSF',QSM_isLambdaCSF,'QSM_lambdaCSF',QSM_lambdaCSF,'eddy',isEddyCorrect,'GPU',isGPU);

        case 'Phase unwrapping'
            % Core of phase unwrapping only 
            UnwrapPhaseMacroIOWrapper(inputDir,outputDir,'FSLBet',isBET,'mask',maskFullName,'unwrap',phaseUnwrap,...
                'Subsampling',subsampling,'exclude_threshold',excludeMaskThreshold,'eddy',isEddyCorrect,'GPU',isGPU);

        case 'Background field removal'
            % core of background field removal only
            BackgroundRemovalMacroIOWrapper(inputDir,outputDir,'mask',maskFullName,...
                'BFR',BFR,'refine',refine,'BFR_tol',BFR_tol,...
                'depth',BFR_depth,'peel',BFR_peel,'BFR_iteration',BFR_iteration,'BFR_padsize',BFR_padSize,...
                'BFR_radius',BFR_radius,'BFR_alpha',BFR_alpha,'BFR_threshold',BFR_threshold,'GPU',isGPU);

        case 'QSM'
            % core of QSM only
            QSMMacroIOWrapper(inputDir,outputDir,'mask',maskFullName,...
                'QSM',QSM_method,'QSM_threshold',QSM_threshold,'QSM_lambda',QSM_lambda,'QSM_optimise',QSM_optimise,...
                'QSM_tol',QSM_tol,'QSM_iteration',QSM_maxiter,'QSM_tol1',QSM_tol1,'QSM_tol2',QSM_tol2,...
                'QSM_padsize',QSM_padsize,'QSM_mu',QSM_mu1,'QSM_mu2',QSM_mu2,QSM_solver,QSM_constraint,...
                'QSM_zeropad',QSM_zeropad,'QSM_wData',QSM_wData,'QSM_wGradient',QSM_wGradient,'QSM_radius',QSM_radius,...
                'QSM_isSMV',QSM_isSMV,'QSM_merit',QSM_merit,'QSM_isLambdaCSF',QSM_isLambdaCSF,'QSM_lambdaCSF',QSM_lambdaCSF,'GPU',isGPU);
    end
catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end

% generate a log file
GenerateLogFile(h.StepsPanel.dataIO.Parent.Title);

% re-enable the pushbutton
set(source,'Enable','on');

function GenerateLogFile(tab)

switch tab
    case 'One-stop QSM processing'
        fid = fopen([outputDir filesep 'qsm_hub.log'],'w');
        fprintf(fid,'QSMHub(''%s'',''%s'',''FSLBet'',%i,''mask'',''%s'',''unwrap'',''%s'',...\n',inputDir,outputDir,isBET,maskFullName,phaseUnwrap);
        fprintf(fid,'''Subsampling'',%i,''BFR'',''%s'',''refine'',%i,''BFR_tol'',%g,...\n',subsampling,BFR,refine,BFR_tol);
        fprintf(fid,'''depth'',%i,''peel'',%i,''BFR_iteration'',%i,''BFR_padsize'',%i,...\n',BFR_depth,BFR_peel,BFR_iteration,BFR_padSize);
        fprintf(fid,'''BFR_radius'',[%s],''BFR_alpha'',%g,''BFR_threshold'',%g,...\n',num2str(BFR_radius),BFR_alpha,BFR_threshold);
        fprintf(fid,'''QSM'',''%s'',''QSM_threshold'',%g,''QSM_lambda'',%g,''QSM_optimise'',%i,...\n',QSM_method,QSM_threshold,QSM_lambda,QSM_optimise);
        fprintf(fid,'''QSM_tol'',%g,''QSM_iteration'',%i,''QSM_tol1'',%g,''QSM_tol2'',%g,...\n',QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2);
        fprintf(fid,'''QSM_padsize'',[%s],''QSM_mu'',%g,''QSM_mu2'',%g,''%s'',''%s'',''exclude_threshold'',%i,...\n',num2str(QSM_padsize),QSM_mu1,QSM_mu2,QSM_solver,QSM_constraint,excludeMaskThreshold);
        fprintf(fid,'''QSM_zeropad'',%i,''QSM_wData'',%g,''QSM_wGradient'',%g,''QSM_radius'',%i,...\n',QSM_zeropad,QSM_wData,QSM_wGradient,QSM_radius);
        fprintf(fid,'''QSM_isSMV'',%i,''QSM_merit'',%i,''QSM_isLambdaCSF'',%g,''QSM_lambdaCSF'',%g,''eddy'',%i,''GPU'',%i);\n',QSM_isSMV,QSM_merit,QSM_isLambdaCSF,QSM_lambdaCSF,isEddyCorrect,isGPU);
        fclose(fid);
        
    case 'Phase unwrapping'
        fid = fopen([outputDir filesep 'qsm_hub.log'],'w');
        fprintf(fid,'UnwrapPhaseMacroIOWrapper(''%s'',''%s'',''FSLBet'',%i,''mask'',''%s'',''unwrap'',''%s'',...\n',inputDir,outputDir,isBET,maskFullName,phaseUnwrap);
        fprintf(fid,'''Subsampling'',%i,''exclude_threshold'',%i,''eddy'',%i,''GPU'',%i);\n',subsampling,excludeMaskThreshold,isEddyCorrect,isGPU);
        fclose(fid);
    
    case 'Background field removal'
        fid = fopen([outputDir filesep 'qsm_hub.log'],'w');
        fprintf(fid,'BackgroundRemovalMacroIOWrapper(''%s'',''%s'',''mask'',''%s'',...\n',inputDir,outputDir,maskFullName);
        fprintf(fid,'''BFR'',''%s'',''refine'',%i,''BFR_tol'',%g,...\n',BFR,refine,BFR_tol);
        fprintf(fid,'''depth'',%i,''peel'',%i,''BFR_iteration'',%i,''BFR_padsize'',%i,...\n',BFR_depth,BFR_peel,BFR_iteration,BFR_padSize);
        fprintf(fid,'''BFR_radius'',[%s],''BFR_alpha'',%g,''BFR_threshold'',%g,''GPU'',%i);\n',num2str(BFR_radius),BFR_alpha,BFR_threshold,isGPU);
        fclose(fid);
    
    case 'QSM'
        fid = fopen([outputDir filesep 'qsm_hub.log'],'w');
        fprintf(fid,'qsmMacroIOWrapper(''%s'',''%s'',''mask'',''%s'',...\n',inputDir,outputDir,maskFullName);
        fprintf(fid,'''QSM'',''%s'',''QSM_threshold'',%g,''QSM_lambda'',%g,''QSM_optimise'',%i,...\n',QSM_method,QSM_threshold,QSM_lambda,QSM_optimise);
        fprintf(fid,'''QSM_tol'',%g,''QSM_iteration'',%i,''QSM_tol1'',%g,''QSM_tol2'',%g,...\n',QSM_tol,QSM_maxiter,QSM_tol1,QSM_tol2);
        fprintf(fid,'''QSM_padsize'',[%s],''QSM_mu'',%g,''QSM_mu2'',%g,''%s'',''%s'',...\n',num2str(QSM_padsize),QSM_mu1,QSM_mu2,QSM_solver,QSM_constraint);
        fprintf(fid,'''QSM_zeropad'',%i,''QSM_wData'',%g,''QSM_wGradient'',%g,''QSM_radius'',%i,...\n',QSM_zeropad,QSM_wData,QSM_wGradient,QSM_radius);
        fprintf(fid,'''QSM_isSMV'',%i,''QSM_merit'',%i,''QSM_isLambdaCSF'',%g,''QSM_lambdaCSF'',%g,''GPU'',%i);\n',QSM_isSMV,QSM_merit,QSM_isLambdaCSF,QSM_lambdaCSF,isGPU);
        fclose(fid);

end

end

end

%% Utility functions

% 1. Get qsm_hub header
function ButtonOpen_Utility_getHeader_Callback(source,eventdata,field)
% get input file/directory for getHeader utility function
global h

switch field
    case 'te'
        % te file can be text file or mat file
        [tefileName,pathDir] = uigetfile({'*.mat;*.txt'},'Select TE file');
        
        % display file directory
        if pathDir ~= 0
            set(h.Utility.getHeader.edit.teFile,'String',fullfile(pathDir,tefileName));
        end

    case 'nitfi'
        % read NIfTI file 
        [nitfiName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select mask file');

        if pathDir ~= 0
            set(h.Utility.getHeader.edit.niftiInput,    'String',fullfile(pathDir,nitfiName));
            % automatically set default output field
            set(h.Utility.getHeader.edit.outputDir,     'String',pathDir);
        end
        
    case 'dicom'
        % get DICOM directory
        pathDir = uigetdir;

        if pathDir ~= 0
            % set input edit field for display
            set(h.Utility.getHeader.edit.dicomInput,    'String',pathDir);
            % automatically set default output field
            set(h.Utility.getHeader.edit.outputDir,     'String',pathDir);
        end
        
    case 'output'
        
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.Utility.getHeader.edit.outputDir,     'String',pathDir);
        end
end

end

function PushbuttonRun_Utility_getHeader_Callback(source,eventdata)
% Callback function to detect and save header functino for qsm_hub

global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

try 

% get DICOM directory (if any)
dicomDir    = get(h.Utility.getHeader.edit.dicomInput,'String'); 
% get output directory, assume the same directory as input directory/file
outputDir = get(h.Utility.getHeader.edit.outputDir, 'String');


% if no dicom directory detected then get te from user input
if isempty(dicomDir)
    
    % get NIfTI file (if any)
    niftiFile   = get(h.Utility.getHeader.edit.niftiInput, 'String');
    
    if isempty(niftiFile)
        % This function requires input file/directory
        error('Please specify a DICOM directory or a NIfTI file.');
    else
        
        % get user input magnetic field strength
        b0          = str2double(get(h.Utility.getHeader.edit.userB0, 'String'));
        % get user input magnetic field direction
        b0dir       = str2num(get(h.Utility.getHeader.edit.userB0dir, 'String'));
        % get user input voxel size
        voxelSize   = str2num(get(h.Utility.getHeader.edit.userVoxelSize, 'String'));

        % check validity of input voxel size
        if length(voxelSize) ~= 3 && ~isempty(voxelSize) && isempty(dicomDir)
            error(['qsm_hub currently works with 3D data only. ' 
                   'Please specify the voxel size in all three dimensions']);
        end

        % check validity of input B0 direction
        if length(b0dir) ~= 3 && ~isempty(b0dir) && isempty(dicomDir)
            error('The B0 direction has 3 elements [x,y,z]');
        end

        % str2double returns NaN if input field is empty, replace it by []
        if isnan(b0)
            b0=[];
        end

        % try to get TE variable in this stage
        teFullName = get(h.Utility.getHeader.edit.teFile, 'String');
        teUserInput = get(h.Utility.getHeader.edit.userTE, 'String');
        if ~isempty(teFullName)
            % get data type
            [~,~,ext] = fileparts(teFullName);

            if strcmpi(ext,'.mat')
                % if mat file then try to load 'TE' directly
                try load(teFullName,'TE');  catch; error('No variable named TE.'); end
            else
                % if text file the try to read the TEs line by line
                TE = readTEfromText(teFullName);
                TE = TE(:);
            end
        else
            % read user input array
            TE = str2num(teUserInput);
            TE = TE(:);
        end
        % in case TE cannot be found
        if isempty(TE) || isnan(TE(1))
            error('Incorrect TE format.');
        end
        
        % specify the input is NIfTI file
        ExportQMSHubHeaderIOWrapper(niftiFile,outputDir,...
            b0,b0dir,voxelSize,TE,'nifti');
    end
else
    
    % specify the input is DICOM directory
    ExportQMSHubHeaderIOWrapper(dicomDir,outputDir,[],[],[],[],'dicom');
    
end

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end

disp('Done!');

% re-enable the pushbutton 
set(source,'Enable','on');

end

% 2. Get lateral ventricle mask 
function ButtonOpen_Utility_csfMask_Callback(source,eventdata,field)
% get input file/directory for getHeader utility function
global h

switch field
    case 'te'
        % te file can be text file or mat file
        [tefileName,pathDir] = uigetfile({'*.mat;*.txt'},'Select TE file');
        
        % display file directory
        if pathDir ~= 0
            set(h.Utility.csfMask.edit.teFile,'String',fullfile(pathDir,tefileName));
        end

    case 'nitfi'
        % read NIfTI file 
        [nitfiName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select NIfTI file');

        if pathDir ~= 0
            set(h.Utility.csfMask.edit.niftiInput,    'String',fullfile(pathDir,nitfiName));
            % automatically set default output field
            set(h.Utility.csfMask.edit.outputDir,     'String',[pathDir filesep 'output']);
        end
        
    case 'mask'
        % read NIfTI file 
        [nitfiName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select mask file');

        if pathDir ~= 0
            set(h.Utility.csfMask.edit.maskInput,    'String',fullfile(pathDir,nitfiName));
        end
        
    case 'output'
        
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.Utility.csfMask.edit.outputDir,     'String',pathDir);
        end
end

end

function PushbuttonRun_Utility_csfMask_Callback(source,event)
% Callback function to get lateral ventricle mask

global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

try

% get NIfTI file (if any)
niftiFile   = get(h.Utility.csfMask.edit.niftiInput, 'String');
% get output directory, assume the same directory as input directory/file
outputDir = get(h.Utility.csfMask.edit.outputDir, 'String');
% mask file
maskFullName    = get(h.Utility.csfMask.edit.maskInput,'String');

if isempty(niftiFile)
    % This function requires input file/directory
    error('Please specify a NIfTI file.');
else

    % try to get TE variable in this stage
    teFullName = get(h.Utility.csfMask.edit.teFile, 'String');
    teUserInput = get(h.Utility.csfMask.edit.userTE, 'String');
    if ~isempty(teFullName)
        % get data type
        [~,~,ext] = fileparts(teFullName);

        if strcmpi(ext,'.mat')
            % if mat file then try to load 'TE' directly
            try load(teFullName,'TE');  catch; error('No variable named TE.'); end
        else
            % if text file the try to read the TEs line by line
            TE = readTEfromText(teFullName);
            TE = TE(:);
        end
    else
        % read user input array
        TE = str2num(teUserInput);
        TE = TE(:);
    end
    % in case TE cannot be found
    if isempty(TE) || isnan(TE(1))
        error('Incorrect TE format.');
    end
    
    GetVentricleMaskIOWrapper(niftiFile,outputDir,maskFullName,TE);
    
end

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end
    

% Disable the pushbutton to prevent doubel click
set(source,'Enable','on');

end

