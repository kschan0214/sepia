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

% function that create the method specific panel
% order must match the order in variable 'methodQSMName'
function_QSM_method_panel = {'sepia_handle_panel_qsm_TKD',...
                             'sepia_handle_panel_qsm_CFS',...
                             'sepia_handle_panel_qsm_NDI',...
                             'sepia_handle_panel_qsm_STIiLSQR',...
                             'sepia_handle_panel_qsm_iLSQR',...
                             'sepia_handle_panel_qsm_FANSI',...
                             'sepia_handle_panel_qsm_Star',...
                             'sepia_handle_panel_qsm_MEDI',...
                             };
                             % in future, add panel of new method here

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
    subwidth(1) = 0.5;
    subwidth(2) = 1-subwidth(1);
    % text|popup pair: select method
    h.qsm.text.qsm = uicontrol('Parent',h.StepsPanel.qsm,'Style','text','String','Method:',...
        'units','normalized','position',[left(1) 0.85 width*subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Select QSM algorithm');
    h.qsm.popup.qsm = uicontrol('Parent',h.StepsPanel.qsm,'Style','popup',...
        'String',methodQSMName,...
        'units','normalized','position',[left(1)+width*subwidth(1) 0.85 width*subwidth(2) height]) ; 
    
    % text|popup pair: select reference tissue
    h.qsm.text.tissue = uicontrol('Parent',h.StepsPanel.qsm,'Style','text','String','Reference tissue:',...
        'units','normalized','position',[left(2) 0.85 width*subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Region used to normalise the magnetic susceptibility map');
    h.qsm.popup.tissue = uicontrol('Parent',h.StepsPanel.qsm,'Style','popup',...
        'String',tissueName,...
        'units','normalized','position',[left(2)+width*subwidth(1) 0.85 width*subwidth(2) height]) ; 
    
%% create control panel

% define position and size of all method panels
position_child = [0.01 0.04 0.95 0.75];

% construct all method panels
for k = 1:length(function_QSM_method_panel)
    h = feval(function_QSM_method_panel{k},h.StepsPanel.qsm,h,position_child);
end

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