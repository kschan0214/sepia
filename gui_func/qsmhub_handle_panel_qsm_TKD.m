%% h = qsmhub_handle_panel_qsm_TKD(hParent,h,position)
%
% Input
% --------------
% hParent       : parent handle of this panel
% h             : global structure contains all handles
% position      : position of this panel
%
% Output
% --------------
% h             : global structure contains all new and other handles
%
% Description: This GUI function creates a panel for LBV method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date last modified: 
%
%
function h = qsmhub_handle_panel_qsm_TKD(hParent,h,position)

%% Parent handle of TKD panel children

h.qsm.panel.TKD = uipanel(hParent,...
    'Title','Thresholded k-space division (TKD)',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','on');

%% Children of TKD panel
    
    % text|edit field pair: threshold
    h.qsm.TKD.text.threshold = uicontrol('Parent',h.qsm.panel.TKD,...
        'Style','text',...
        'String','Threshold (0-1):',...
        'units','normalized','position',[0.01 0.75 0.23 0.2],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','K-space threshold');
    h.qsm.TKD.edit.threshold = uicontrol('Parent',h.qsm.panel.TKD,...
        'Style','edit',...
        'String','0.15',...
        'units','normalized','position',[0.25 0.75 0.2 0.2],...
        'backgroundcolor','white');

end