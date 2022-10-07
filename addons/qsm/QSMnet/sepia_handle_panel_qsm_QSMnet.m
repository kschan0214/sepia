%% h = sepia_handle_panel_qsm_QSMnet(hParent,h,position)
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
% Description: This GUI function creates a panel for NDI method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 October 2022
% Date modified: 
%
%
function h = sepia_handle_panel_qsm_QSMnet(hParent,h,position)

%% set default values

%% Tooltips

%% layout of the panel

%% Parent handle of QSMnet panel children

h.qsm.panel.QSMnet = uipanel(hParent,...
    'Title','QSMnet+',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of QSMnet panel

end