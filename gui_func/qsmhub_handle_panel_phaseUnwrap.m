%% function output = function_name(input)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
%
function h = qsmhub_handle_panel_phaseUnwrap(hParent,hFig,h,position)
h.StepsPanel.phaseUnwrap = uipanel(hParent,'Title','Total field recovery and phase unwrapping',...
    'position',[position(1) position(2) 0.95 0.2]);
    h.phaseUnwrap.text.phaseUnwrap = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,'Style','text','String','Phase unwrapping:',...
        'units','normalized','position',[0.01 0.84 0.3 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(hFig,'color'),...
        'tooltip','Select phase unwrapping algorithm');
    h.phaseUnwrap.popup.phaseUnwrap = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,'Style','popup',...
        'String',{'Laplacian','Laplacian STI suite','Jena','Region growing','Graphcut'},...
        'units','normalized','position',[0.31 0.84 0.4 0.15]); 
%     h.phaseUnwrap.text.unit = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,'Style','text','String','Output unit:',...
%         'units','normalized','position',[0.01 0.65 0.3 0.15],...
%         'HorizontalAlignment','left',...
%         'backgroundcolor',get(hFig,'color'),...
%         'tooltip','Select unwrapped phase map unit');
%     h.phaseUnwrap.popup.unit = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,'Style','popup',...
%         'String',{'radHz','ppm','rad','Hz'},...
%         'units','normalized','position',[0.31 0.65 0.4 0.15]);
    h.phaseUnwrap.checkbox.eddyCorrect = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,'Style','checkbox','String','Bipolar readout eddy current correction',...
        'units','normalized','Position',[0.01 0.41 1 0.2],...
        'backgroundcolor',get(hFig,'color'),...
        'tooltip','Correct the inconsistent phase between odd and even echoes');
    h.phaseUnwrap.checkbox.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,'Style','checkbox','String','Exclude unreliable voxels',...
        'units','normalized','position',[0.01 0.20 0.6 0.2],...
        'backgroundcolor',get(hFig,'color'));
    h.phaseUnwrap.text.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,'Style','text','String','Threshold (0-1):',...
        'units','normalized','position',[0.01 0.02 0.3 0.15],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(hFig,'color'),...
        'tooltip','Threshold to exclude unreliable voxels [0,1]. Smaller value gives smaller brain volume.');
    h.phaseUnwrap.edit.excludeMask = uicontrol('Parent',h.StepsPanel.phaseUnwrap ,'Style','edit',...
            'String','0.0009',...
            'units','normalized','position',[0.31 0.02 0.3 0.15],...
            'backgroundcolor','white',...
            'Enable','off');
end