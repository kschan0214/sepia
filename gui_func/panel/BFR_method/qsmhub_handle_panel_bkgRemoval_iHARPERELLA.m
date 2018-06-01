%% h = qsmhub_handle_panel_bkgRemoval_iHARPERELLA(hParent,h,position)
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
% Description: This GUI function creates a panel for iHARPERELLA 
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date last modified: 
%
%
function h = qsmhub_handle_panel_bkgRemoval_iHARPERELLA(hParent,h,position)

%% Parent handle of iHARPERELLA panel

h.bkgRemoval.panel.iHARPERELLA = uipanel(hParent,...
    'Title','iHARPERELLA',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Chrild of iHARPERELLA panel

    % text|edit field pair: maximum number of iterations
    h.bkgRemoval.iHARPERELLA.text.maxIter = uicontrol('Parent',h.bkgRemoval.panel.iHARPERELLA,...
        'Style','text',...
        'String','Max. iterations:',...
        'units','normalized','position',[0.01 0.75 0.2 0.2],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Maximum iteration allowed');
    h.bkgRemoval.iHARPERELLA.edit.maxIter = uicontrol('Parent',h.bkgRemoval.panel.iHARPERELLA,...
        'Style','edit',...
        'String','100',...
        'units','normalized','position',[0.25 0.75 0.2 0.2],...
        'backgroundcolor','white');

%% set callbacks
set(h.bkgRemoval.iHARPERELLA.edit.maxIter,	'Callback', {@EditInputMinMax_Callback,1,0});

end