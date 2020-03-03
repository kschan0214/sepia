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
% Date modified: 3 March 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_qsm(hParent,h,position)
% set up method name displayed on GUI
methodName = {'TKD','Closed-form solution','NDI','STI suite iLSQR','iLSQR','FANSI','Star-QSM','MEDI'};
tissueName = {'None','Brain mask','CSF'};

%% layout of the panel
% define maximum level of options and spacing between options
ncol    = 2; % 2 columns in the panel
rspacing = 0.01;
width   = (1-(ncol+1)*rspacing)/ncol;
left    = (rspacing:width+rspacing:(width+rspacing)*ncol);

% Set parent of qsm panel
h.StepsPanel.qsm = uipanel(hParent,...
    'Title','QSM','backgroundcolor',get(h.fig,'color'),...
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
        'String',methodName,...
        'units','normalized','position',[left(1)+width*subwidth(1) 0.85 width*subwidth(2) height]) ; 
    
    % text|popup pair: select reference tissue
    h.qsm.text.tissue = uicontrol('Parent',h.StepsPanel.qsm,'Style','text','String','Reference tissue:',...
        'units','normalized','position',[left(2) 0.85 width*subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Region used to normalised the magnetic susceptibility map');
    h.qsm.popup.tissue = uicontrol('Parent',h.StepsPanel.qsm,'Style','popup',...
        'String',tissueName,...
        'units','normalized','position',[left(2)+width*subwidth(1) 0.85 width*subwidth(2) height]) ; 
    
%% create control panel

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
    
    % NDI
    h = sepia_handle_panel_qsm_NDI(h.StepsPanel.qsm,h,position_child);
    
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
        
    case 'NDI'
        set(h.qsm.panel.NDI,        'Visible','on');

    % in the future, add new method here

end

end