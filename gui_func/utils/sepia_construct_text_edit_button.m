%% [h_text,h_edit] = sepia_construct_text_edit(parent,fieldString,defaultValue,pos,subratio)
%
% Input
% --------------
% parent        : parent of the constructs
% fieldString   : text to be displayed in the 'text' field
% defaultValue  : default value in the 'edit' field
% buttonDisp    : string or image displayed in the 'push button' field
% pos           : position of the constructs, [left,bottom,total_width,height]
% wratio        : normalised width proportion of all fields
%
% Output
% --------------
% h_text        : handle of the text field
% h_edit        : handle of the edit field
% h_putton      : handle of the push button field
%
% Description: constructs a 'text|edit|button' pair on the GUIs
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 5 April 2020
% Date modified:
%
%
function [h_text,h_edit,h_putton] = sepia_construct_text_edit_button(parent,fieldString,defaultValue,buttonDisp,pos,wratio)

if nargin < 6
    wratio = [0.3,0.65,0.05];
end
wratio = wratio ./ sum(wratio);

width = pos(3) * wratio;

h_text = uicontrol('Parent',parent ,...
    'Style','text',...
    'String',fieldString,...
    'units','normalized','position',[pos(1) pos(2) width(1) pos(4)],...
    'HorizontalAlignment','left',...
    'backgroundcolor',get(gcf,'color'));

h_edit = uicontrol('Parent',parent ,...
    'Style','edit',...
    'String',num2str(defaultValue),...
    'units','normalized','position',[pos(1)+width(1) pos(2) width(2) pos(4)],...
    'backgroundcolor','white');

h_putton = uicontrol('Parent',parent,...
    'Style','pushbutton',...
    'units','normalized','position',[pos(1)+sum(width(1:2)) pos(2) width(3) pos(4)],...
    'backgroundcolor','white');    

if ischar(buttonDisp)
    set(h_putton, 'String', buttonDisp);
else
    set(h_putton, 'CData', buttonDisp);
end
    
end