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
defaultMFGThreshold  = 0.5;

tooltip.QSM.panel.method    = 'Select a QSM algorithm'; 
tooltip.QSM.panel.reference	= 'Region used to normalise the magnetic susceptibility map';
tooltip.QSM.panel.twopass	= 'Two pass masking mask refinement strategy to use';

%% layout of the panel
% define maximum level of options and spacing between options
ncol    = 2; % 2 columns in the panel
rspacing = 0.01;
width   = (1-(ncol+1)*rspacing)/ncol;
left    = (rspacing:width+rspacing:(width+rspacing)*ncol);

% Set parent of qsm panel
h.StepsPanel.qsm = uipanel(hParent,...
    'Title','QSM','backgroundcolor',get(h.fig,'color'),...
    'fontweight', 'bold',...
    'position',[position(1) position(2) 0.95 0.25]);

%% design of this panel
    
    height = 0.12;
    wratio = 0.5;

    % col 1
    % text|popup pair: select method
    [h.qsm.text.qsm,h.qsm.popup.qsm] = sepia_construct_text_popup(...
        h.StepsPanel.qsm,'Method:', methodQSMName, [left(1) 0.85 width height], wratio);
    
    % utility function related to two pass masking
    h.qsm.text.twopass = uicontrol('Parent',h.StepsPanel.qsm,...
        'Style','text',...
        'String','Perform two pass masking',...
        'units','normalized','position',[left(1) 0.03 width*0.3 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.qsm.popup.twopass = uicontrol('Parent',h.StepsPanel.qsm ,...
        'Style','popup',...
        'String',methodTwoPassName,...
        'Enable','on',...
        'units','normalized','position',[left(1)+width*0.3 0.03 width*0.35 height]); 
    h.qsm.text.lambda = uicontrol('Parent',h.StepsPanel.qsm ,...
        'Style','text','String','λ:',...
        'units','normalized','position',[left(1)+width*0.7 0.03 width*0.1 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.qsm.edit.lambda = uicontrol('Parent',h.StepsPanel.qsm ,...
        'Style','edit',...
        'String',num2str(defaultMFGThreshold),...
        'units','normalized','position',[left(1)+width*0.75 0.03 width*0.1 height],...
        'backgroundcolor','white',...
        'Enable','on');
    h.qsm.slider.lambda = uicontrol('Parent',h.StepsPanel.qsm ,...
        'Style','slider',...
        'Value',defaultMFGThreshold,...
        'units','normalized','position',[left(1)+width*0.85 0.03 0.01 height],...
        'Max',4, 'Min',0,'SliderStep',[0.25 0.50],...
        'Enable','on');

    % col 2
    % text|popup pair: select reference tissue
    [h.qsm.text.tissue,h.qsm.popup.tissue] = sepia_construct_text_popup(...
        h.StepsPanel.qsm,'Reference tissue:', tissueName, [left(2) 0.85 width height], wratio);

    
%% create control panel

% define position and size of all method panels
position_child = [0.01 0.2 0.95 0.65];

% construct all method panels
for k = 1:length(function_QSM_method_panel)
    h = feval(function_QSM_method_panel{k},h.StepsPanel.qsm,h,position_child);
end

%% set tooltip
set(h.qsm.text.qsm,     'Tooltip',tooltip.QSM.panel.method);
set(h.qsm.text.tissue,  'Tooltip',tooltip.QSM.panel.reference);
set(h.qsm.text.twopass, 'Tooltip',tooltip.QSM.panel.twopass);

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

% switch off all panels
fields = fieldnames(h.qsm.panel);
for kf = 1:length(fields)
    set(h.qsm.panel.(fields{kf}),   'Visible','off');
end

% switch on only target panel
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
    
    case methodTwoPassName{1}
    % get slider value and update the edit field
    set(h.qsm.edit.lambda,    'String', 0.5);
    set(h.qsm.slider.lambda,  'Value',  0.5);
    set(h.qsm.edit.lambda,    'enable', 'on');
    set(h.qsm.slider.lambda,  'enable', 'on');

    case methodTwoPassName{2}
    % get slider value and update the edit field
    set(h.qsm.edit.lambda,    'enable', 'off');
    set(h.qsm.slider.lambda,  'enable', 'off');
end

end
