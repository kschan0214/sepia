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
% Date last modified: 1 June 2018
%
%
function h = sepia_handle_panel_qsm(hParent,h,position)
% set up method name displayed on GUI
methodName = {'TKD','Closed-form solution','STI suite iLSQR','iLSQR','FANSI','Star-QSM','MEDI'};

% Set parent of qsm panel
h.StepsPanel.qsm = uipanel(hParent,...
    'Title','QSM','backgroundcolor',get(h.fig,'color'),...
    'position',[position(1) position(2) 0.95 0.25]);

%% design of this panel

    % text|popup pair: select method
    h.qsm.text.qsm = uicontrol('Parent',h.StepsPanel.qsm,'Style','text','String','Method:',...
        'units','normalized','position',[0.01 0.85 0.15 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Select QSM algorithm');
    h.qsm.popup.qsm = uicontrol('Parent',h.StepsPanel.qsm,'Style','popup',...
        'String',methodName,...
        'units','normalized','position',[0.31 0.85 0.4 0.1]) ; 
    
% create control panel

% define position and size of all method panels
position_child = [0.01 0.04 0.95 0.75];

    % TKD    
    h = sepia_handle_panel_qsm_TKD(h.StepsPanel.qsm,h,position_child);

    % Closed-form solution  
    h = sepia_handle_panel_qsm_CFS(h.StepsPanel.qsm,h,position_child);

    % iLSQR  
    h = sepia_handle_panel_qsm_iLSQR(h.StepsPanel.qsm,h,position_child);
        
    % STI suite iLSQR  
    h = sepia_handle_panel_qsm_STIiLSQR(h.StepsPanel.qsm,h,position_child);
    
    % FANSI
    h = sepia_handle_panel_qsm_FANSI(h.StepsPanel.qsm,h,position_child);
        
    % Star-QSM
    h = sepia_handle_panel_qsm_Star(h.StepsPanel.qsm,h,position_child);

    % MEDI
    h = sepia_handle_panel_qsm_MEDI(h.StepsPanel.qsm,h,position_child);
    
    % in future, add panel of new method here

% set callback
set(h.qsm.popup.qsm, 'Callback', {@PopupQSM_Callback,h});

end

%% Callback function
function PopupQSM_Callback(source,eventdata,h)
% display corresponding QSM method's panel

% global h

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

    case 'Star-QSM'
        set(h.qsm.panel.Star,       'Visible','on');

    case 'MEDI'
        set(h.qsm.panel.MEDI,       'Visible','on');

    % in the future, add new method here

end

end