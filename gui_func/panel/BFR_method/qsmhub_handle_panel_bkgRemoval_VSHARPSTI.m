%% h = qsmhub_handle_panel_bkgRemoval_VSHARPSTI(hParent,h,position)
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
% Description: This GUI function creates a panel for STI suite VSHARP 
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date last modified: 
%
%
function h = qsmhub_handle_panel_bkgRemoval_VSHARPSTI(hParent,h,position)

defaultSMVSize = 12;

%% Parent handle of VSHARPSTI panel children

h.bkgRemoval.panel.VSHARPSTI = uipanel(hParent,...
    'Title','STI suite Variable SHARP',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Child of VSHARPSTI panel

    % text|edit field pair: SMV size
    h.bkgRemoval.VSHARPSTI.text.smvSize = uicontrol('Parent',h.bkgRemoval.panel.VSHARPSTI,...
        'Style','text',...
        'String','SMV size (mm):',...
        'units','normalized','position',[0.01 0.75 0.2 0.2],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.bkgRemoval.VSHARPSTI.edit.smvSize = uicontrol('Parent',h.bkgRemoval.panel.VSHARPSTI,...
        'Style','edit',...
        'String',num2str(defaultSMVSize),...
        'units','normalized','position',[0.25 0.75 0.2 0.2],...
        'backgroundcolor','white');

%% set callbacks
set(h.bkgRemoval.VSHARPSTI.edit.smvSize,    'Callback', {@EditInputMinMax_Callback,defaultSMVSize,1,0});

end