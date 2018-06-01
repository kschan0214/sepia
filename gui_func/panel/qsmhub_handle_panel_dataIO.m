%% h = qsmhub_handle_panel_dataIO(hParent,h,position)
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
% Description: This GUI function creates a panel for I/O control
%
% Kwok-Shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 April 2018
% Date last modified: 30 May 2018
%
%
function h = qsmhub_handle_panel_dataIO(hParent,h,position)

% define maximum level of options and spacing between options
nlevel = 5;
spacing = 0.02;
height = (1-(nlevel+1)*spacing)/nlevel;
button = (height+spacing:height+spacing:(height+spacing)*nlevel) - height;

% Parent of dataIO panel
h.StepsPanel.dataIO = uipanel(hParent,'Title','I/O',...
    'Position',[position(1) position(2) 0.95 0.2],...
    'backgroundcolor',get(h.fig,'color'));

    % input directory
    h.dataIO.text.input = uicontrol(h.StepsPanel.dataIO,...
        'Style','text','String','Input directory:',...
        'Units','normalized','Position', [0.01 button(5) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Input directory contains DICOM (both magnitude and phase files under the same directory) or NIfTI (*phase*.nii* and *magn*.nii*) files');
    h.dataIO.edit.input = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(5) 0.68 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.dataIO.button.input = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 button(5) 0.1 height],...
        'backgroundcolor','white');

    % output basename
    h.dataIO.text.output = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','text','String','Output basename:',...
        'units','normalized','Position',[0.01 button(4) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Output basename. Default output basename is the input directory followed by /output/squirrel');
    h.dataIO.edit.output = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(4) 0.68 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.dataIO.button.output = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 button(4) 0.1 height],...
        'backgroundcolor','white');

    % invert pahse
    h.dataIO.checkbox.invertPhase = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','checkbox','String','Invert phase data',...
        'units','normalized','Position',[0.01 button(3) 0.5 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Reverse the direction of frequency shift. This could be useful if your previous QSM results show negative values in deep gray matter structure.');

    % brain extraction
    h.dataIO.checkbox.brainExtraction = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','checkbox','String','FSL brain extraction',...
        'units','normalized','Position',[0.01 button(2) 0.5 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Use FSL bet for brain extraction');

    % brain mask
    h.dataIO.text.maskdir = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','text','String','Brain mask file:',...
        'units','normalized','Position',[0.01 button(1) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Specify mask file (NIfTI format). If the input directory contains NIfTI file with *mask* in the filename then it will be read automatically');
    h.dataIO.edit.maskdir = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(1) 0.68 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.dataIO.button.maskdir = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 button(1) 0.1 height],...
        'backgroundcolor','white');

%% set callback function
% checkbox/edit pair
set(h.dataIO.checkbox.brainExtraction,  'Callback', {@CheckboxBrainExtraction_Callback,h});
% open/edit pari
set(h.dataIO.button.input,           	'Callback', {@ButtonOpen_Callback,h,'input'});
set(h.dataIO.button.output,         	'Callback', {@ButtonOpen_Callback,h,'output'});
set(h.dataIO.button.maskdir,          	'Callback', {@ButtonOpen_Callback,h,'mask'});
    
end