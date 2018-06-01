%% h = qsmhub_handle_panel_qsm_Star(hParent,h,position)
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
% Description: This GUI function creates a panel for FANSI method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date last modified: 
%
%
function h = qsmhub_handle_panel_qsm_Star(hParent,h,position)

%% Parent handle of CFS panel children

h.qsm.panel.Star = uipanel(hParent,...
    'Title','Star-QSM',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of CFS panel
    
    % text|edit field pair: pad size
    h.qsm.Star.text.padSize = uicontrol('Parent',h.qsm.panel.Star,...
        'Style','text',...
        'String','Pad size:',...
        'units','normalized','position',[0.01 0.75 0.2 0.2],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.qsm.Star.edit.padSize = uicontrol('Parent',h.qsm.panel.Star,...
        'Style','edit',...
        'String','12',...
        'units','normalized','position',[0.25 0.75 0.2 0.2],...
        'backgroundcolor','white');

%% set callbacks
set(h.qsm.Star.edit.padSize, 'Callback', {@EditInputMinMax_Callback,1,0});

end