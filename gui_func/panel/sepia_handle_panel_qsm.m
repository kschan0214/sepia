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

tooltip.QSM.panel.method    = 'Select a QSM algorithm'; 
tooltip.QSM.panel.reference	= 'Region used to normalise the magnetic susceptibility map';

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
    
    height = 0.1;
    wratio = 0.5;

    % col 1
    % text|popup pair: select method
    [h.qsm.text.qsm,h.qsm.popup.qsm] = sepia_construct_text_popup(...
        h.StepsPanel.qsm,'Method:', methodQSMName, [left(1) 0.85 width height], wratio);
    
    % col 2
    % text|popup pair: select reference tissue
    [h.qsm.text.tissue,h.qsm.popup.tissue] = sepia_construct_text_popup(...
        h.StepsPanel.qsm,'Reference tissue:', tissueName, [left(2) 0.85 width height], wratio);

    
%% create control panel

% define position and size of all method panels
position_child = [0.01 0.04 0.95 0.75];

% construct all method panels
for k = 1:length(function_QSM_method_panel)
    h = feval(function_QSM_method_panel{k},h.StepsPanel.qsm,h,position_child);
end

%% set tooltip
set(h.qsm.text.qsm,     'Tooltip',tooltip.QSM.panel.method);
set(h.qsm.text.tissue,  'Tooltip',tooltip.QSM.panel.reference);

%% set callback
set(h.qsm.popup.qsm, 'Callback', {@PopupQSM_Callback,h});

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