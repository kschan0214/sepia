%% h = qsmhub_handle_panel_bkgRemoval(hParent,h,position)
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
% Description: This GUI function creates a panel for background field
% removal method control
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 April 2018
% Date last modified: 1 June 2018
%
%
function h = qsmhub_handle_panel_bkgRemoval(hParent,h,position)
% set up method name displayed on GUI
methodName = {'LBV','PDF','RESHARP','SHARP','VSHARP STI suite','VSHARP','iHARPERELLA'};

% Set parent of background removal panel
h.StepsPanel.bkgRemoval = uipanel(hParent,...
    'Title','Background field removal',...
    'position',[position(1) position(2) 0.95 0.25],...
    'backgroundcolor',get(h.fig,'color'));

%% design of this panel

    % text|popup pair: select method
    h.bkgRemoval.text.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','text','String','Method:',...
        'units','normalized','position',[0.01 0.85 0.3 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Select background field removal method');
    h.bkgRemoval.popup.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','popup',...
        'String',methodName,...
        'units','normalized','position',[0.31 0.85 0.4 0.1]) ;   
    
    % utility function related to background field removal
    h.bkgRemoval.checkbox.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','checkbox',...
        'String','Remove B1 phase residual by 3D polynomial fit',...
        'units','normalized','position',[0.01 0.01 0.98 0.1],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Roemer coil combination fails to remove B1+ phase and includes some contribution from B1- of the body coil. These can be removed by fitting a 4-th order 3D polynomial on the computed local field map',...
        'value',1);
    
%% create control panel

% define position and size of all method panels
position_child = [0.01 0.15 0.95 0.65];

    % LBV
    h = qsmhub_handle_panel_bkgRemoval_LBV(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
        
    % PDF
    h = qsmhub_handle_panel_bkgRemoval_PDF(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
        
    % SHARP
    h = qsmhub_handle_panel_bkgRemoval_SHARP(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
    
    % RESHARP    
    h = qsmhub_handle_panel_bkgRemoval_RESHARP(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
        
    % VSHARPSTI    
    h = qsmhub_handle_panel_bkgRemoval_VSHARPSTI(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
    
    % VSHARP
    h = qsmhub_handle_panel_bkgRemoval_VSHARP(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);

    % iHARPERELLA    
    h = qsmhub_handle_panel_bkgRemoval_iHARPERELLA(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);

    % in future, add panel of new method here

%% set callback function
set(h.bkgRemoval.popup.bkgRemoval, 'Callback', {@PopupBkgRemoval_Callback,h});
end