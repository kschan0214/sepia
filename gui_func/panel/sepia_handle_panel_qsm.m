%% h = sepia_handle_panel_qsm(hParent,h,position)
%
% Input
% --------------
% hParent       : parent handle of this panel
% hFig          : handle of the GUI
% h             : global structure contains all handles
% position      : position of this panel
%
% Output
% --------------
% h             : global structure contains all new and other handles
%
% Description: This GUI function creates a panel for QSM method control
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 April 2018
% Date modified: 1 June 2018
% Date modified: 5 Juen 2019
% Date modified: 6 March 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_qsm(hParent,h,position)
% set up method name displayed on GUI
sepia_universal_variables;

% Default value
defaultMFGThreshold  = 0.7;

% Default QSM dipole inversion method: prefer FANSI, then MEDI, then the
% LSQR+HEIDI addon (Linux only, bundled at
% <SEPIA_HOME>/../external/HEIDI_SEPIAready), falling back to TKD
% (self-contained) if none of these are available.
toolboxPaths = get_sepia_toolbox_home_paths();
isHEIDIavailable = isunix && ~ismac && exist(fullfile(SEPIA_HOME,'..','external','HEIDI_SEPIAready'),'dir') == 7;
if exist(toolboxPaths.FANSI_HOME,'dir') == 7
    defaultQSMMethod = 'FANSI';
elseif exist(toolboxPaths.MEDI_HOME,'dir') == 7
    defaultQSMMethod = 'MEDI';
elseif isHEIDIavailable
    defaultQSMMethod = 'LSQR+HEIDI';
else
    defaultQSMMethod = 'TKD';
end
defaultQSMIdx = find(strcmpi(methodQSMName, defaultQSMMethod), 1);
if isempty(defaultQSMIdx)
    defaultQSMIdx = 1;
end

tooltip.QSM.panel.method    = 'Select a QSM algorithm';
tooltip.QSM.panel.reference	= 'Region used to normalise the magnetic susceptibility map';
tooltip.QSM.panel.twopass	= 'Two pass masking mask refinement strategy to use';
tooltip.QSM.panel.isHeidi	= 'Reduce streaking artefact using HEIDI. Adjust HEIDI''s parameters in the LSQR+HEIDI panel and switch back to the target dipole inversion method.';

%% layout of the panel
% define maximum level of options and spacing between options
ncol    = 4; % 2 columns in the panel
rspacing = 0.01;
width   = (1-(ncol+1)*rspacing)/ncol;
left    = (rspacing:width+rspacing:(width+rspacing)*ncol);

% Set parent of qsm panel
h.StepsPanel.qsm = uipanel(hParent,...
    'Title','QSM','backgroundcolor',get(h.fig,'color'),...
    'fontweight', 'bold',...
    'position',[position(1) position(2) 0.95 0.25]);

%% design of this panel
    
    height = 0.1;
    wratio = 0.4;

    % col 1
    % text|popup pair: select method
    [h.qsm.text.qsm,h.qsm.popup.qsm] = sepia_construct_text_popup(...
        h.StepsPanel.qsm,'Method:', methodQSMName, [left(1) 0.85 width height], wratio);
    set(h.qsm.popup.qsm, 'Value', defaultQSMIdx);
    
    % col 2
    % utility function related to two pass masking
    ratio1 = 0.35; ratio2 = 0.45; ratio3 = 0.05; ratio4 = 0.1;
    h.qsm.text.twopass = uicontrol('Parent',h.StepsPanel.qsm,...
        'Style','text',...
        'String','2-pass masking:',...
        'units','normalized','position',[left(2) 0.85 width*ratio1 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.qsm.popup.twopass = uicontrol('Parent',h.StepsPanel.qsm ,...
        'Style','popup',...
        'String',methodTwoPassName,...
        'Enable','on',...
        'units','normalized','position',[left(2)+width*ratio1 0.85 width*ratio2 height]); 
    h.qsm.text.lambda = uicontrol('Parent',h.StepsPanel.qsm ,...
        'Style','text','String','λ:',...
        'units','normalized','position',[left(2)+width*(ratio2+ratio1) 0.85 width*ratio3 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.qsm.edit.lambda = uicontrol('Parent',h.StepsPanel.qsm ,...
        'Style','edit',...
        'String',num2str(defaultMFGThreshold),...
        'units','normalized','position',[left(2)+width*(ratio2+ratio1+ratio3) 0.85 width*ratio4 height],...
        'backgroundcolor',sepia_theme_bg_color(),'foregroundcolor',sepia_theme_fg_color(),...
        'Enable','off');
    h.qsm.slider.lambda = uicontrol('Parent',h.StepsPanel.qsm ,...
        'Style','slider',...
        'Value',defaultMFGThreshold,...
        'units','normalized','position',[left(2)+width*(ratio2+ratio1+ratio3+ratio4) 0.85 0.01 height],...
        'Max',4, 'Min',0,'SliderStep',[0.25 0.50],...
        'Enable','off');

    % col 3
    % check box | Streaking artefact reduction using HEIDI
    h.qsm.checkbox.isHeidi = uicontrol('Parent',h.StepsPanel.qsm,...
        'Style','checkbox','String','Streaking reduction by HEIDI',...
        'units','normalized','Position',[left(3) 0.85 width height],...
        'backgroundcolor',get(h.fig,'color'));

    % col 4;
    % text|popup pair: select reference tissue
    [h.qsm.text.tissue,h.qsm.popup.tissue] = sepia_construct_text_popup(...
        h.StepsPanel.qsm,'Reference:', tissueName, [left(4) 0.85 width height], wratio);

    
%% create control panel

% define position and size of all method panels
position_child = [0.01 0.05 0.95 0.75];

% construct all method panels
for k = 1:length(function_QSM_method_panel)
    h = feval(function_QSM_method_panel{k},h.StepsPanel.qsm,h,position_child);
end

% each method panel is constructed with its own hardcoded 'Visible'
% state; sync it here so the one matching the (possibly toolbox-dependent)
% default selection above is the one actually shown
sync_qsm_panel_visibility(h, methodQSMName, defaultQSMMethod);

%% set tooltip
set(h.qsm.text.qsm,     'Tooltip',tooltip.QSM.panel.method);
set(h.qsm.text.tissue,  'Tooltip',tooltip.QSM.panel.reference);
set(h.qsm.text.twopass, 'Tooltip',tooltip.QSM.panel.twopass);
set(h.qsm.checkbox.isHeidi, 'Tooltip',tooltip.QSM.panel.isHeidi);

%% set callback
set(h.qsm.popup.qsm,     'Callback', {@PopupQSM_Callback,h});
set(h.qsm.edit.lambda,   'Callback', {@EditInputMinMax_Callback,defaultMFGThreshold,1,0,10});
set(h.qsm.slider.lambda, 'Callback', {@SliderMGFlambda_Callback,h});
set(h.qsm.popup.twopass, 'Callback', {@PopupTwoPassMasking_Callback,h});

end

%% Callback function
% display corresponding QSM method's panel
function PopupQSM_Callback(source,eventdata,h)

% needs 'methodQSMName' here
sepia_universal_variables;

% get selected QSM method
method = source.String{source.Value,1} ;

% LSQR+HEIDI ships Linux-only compiled binaries (see the platform check in
% Wrapper_QSM_LSQRandHEIDI.m); warn immediately on selection rather than
% only failing once the user has already configured and run the pipeline
if strcmpi(method,'LSQR+HEIDI') && (~isunix || ismac)
    warndlg('LSQR+HEIDI is only supported on Linux systems. Running the pipeline with this method selected on the current platform will fail.', ...
        'Unsupported platform');
end

sync_qsm_panel_visibility(h, methodQSMName, method);

end

% switch on the method panel matching 'method', switch off all others
function sync_qsm_panel_visibility(h, methodQSMName, method)

fields = fieldnames(h.qsm.panel);
for kf = 1:length(fields)
    set(h.qsm.panel.(fields{kf}),   'Visible','off');
end

for k = 1:length(methodQSMName)
    if strcmpi(method,methodQSMName{k})
        set(h.qsm.panel.(fields{k}), 'Visible','on');
        break
    end
end

end

% callback for refine slider
function SliderMGFlambda_Callback(source,eventdata,h)

% get slider value and update the edit field
set(h.qsm.edit.lambda, 'String', num2str(source.Value))

end

% callback for refine edit
function PopupTwoPassMasking_Callback(source,eventdata,h)
sepia_universal_variables;

switch source.String{source.Value}
    
    case methodTwoPassName{3} % for MGF only
        % get slider value and update the edit field
        set(h.qsm.edit.lambda,    'String', num2str(h.qsm.edit.lambda.String)); % KC: bug fix for reading the wrong text
        set(h.qsm.slider.lambda,  'Value',  source.Value);
        set(h.qsm.edit.lambda,    'enable', 'on');
        set(h.qsm.slider.lambda,  'enable', 'on');

    % case methodTwoPassName{2}
    otherwise
        % get slider value and update the edit field
        set(h.qsm.edit.lambda,    'enable', 'off');
        set(h.qsm.slider.lambda,  'enable', 'off');
end

end
