%% [h_text,h_popup] = sepia_construct_text_popup(parent,fieldString,popupMenu,pos,wratio)
%
% Input
% --------------
% parent        : parent of the constructs
% fieldString   : text to be displayed in the 'text' field
% popupMenu     : selection menu in the 'popup' field, cell variable
% pos           : position of the constructs, [left,bottom,total_width,height]
% wratio        : normalised width proportion of the 'text' field, [0,1]
%
% Output
% --------------
% h_text        : handle of the text field
% h_popup       : handle of the popup field
%
% Description: constructs a 'text|popup' pair on the GUIs
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 April 2020
% Date modified:
%
%
function [h_text,h_popup] = sepia_construct_text_popup(parent,fieldString,popupMenu,pos,wratio)

if nargin < 5
    wratio = 0.5;
end

width = pos(3) * [wratio,1-wratio];

h_text = uicontrol('Parent',parent ,...
    'Style','text',...
    'String',fieldString,...
    'units','normalized','position',[pos(1) pos(2) width(1) pos(4)],...
    'HorizontalAlignment','left',...
    'backgroundcolor',get(gcf,'color'));

h_popup = uicontrol('Parent',parent ,...
    'Style','popup',...
    'String',popupMenu,...
    'units','normalized','position',[pos(1)+width(1) pos(2) width(2) pos(4)],...
    'backgroundcolor','white');    
    
end