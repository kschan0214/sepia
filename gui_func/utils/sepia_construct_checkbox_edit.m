%% [h_checkbox,h_edit] = sepia_construct_checkbox_edit(parent,fieldString,defaultValue,pos,wratio)
%
% Input
% --------------
% parent        : parent of the constructs
% fieldString   : text to be displayed in the 'checkbox' field
% defaultValue  : default value in the 'edit' field
% pos           : position of the constructs, [left,bottom,total_width,height]
% wratio        : normalised width proportion of the 'checkbox' field, [0,1]
%
% Output
% --------------
% h_checkbox    : handle of the checkbox field
% h_edit        : handle of the edit field
%
% Description: constructs a 'checkbox|edit' pair on the GUIs, by default
% the 'edit' field is disable
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 April 2020
% Date modified:
%
%
function [h_checkbox,h_edit] = sepia_construct_checkbox_edit(parent,fieldString,defaultValue,pos,wratio)

if nargin < 5
    wratio = 0.5;
end

width = pos(3) * [wratio,1-wratio];

h_checkbox = uicontrol('Parent',parent ,...
    'Style','checkbox',...
    'String',fieldString,...
    'units','normalized','position',[pos(1) pos(2) width(1) pos(4)],...
    'HorizontalAlignment','left',...
    'backgroundcolor',get(gcf,'color'));

h_edit = uicontrol('Parent',parent ,...
    'Style','edit',...
    'String',num2str(defaultValue),...
    'enable','off',...
    'units','normalized','position',[pos(1)+width(1) pos(2) width(2) pos(4)],...
    'backgroundcolor','white');
    
end