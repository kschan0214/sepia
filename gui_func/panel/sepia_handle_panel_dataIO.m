%% h = sepia_handle_panel_dataIO(hParent,h,position)
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
function h = sepia_handle_panel_dataIO(hParent,h,position)

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
    
    % input data1
    h.dataIO.text.inputData1 = uicontrol(h.StepsPanel.dataIO,...
        'Style','text','String','or Phase data:',...
        'Units','normalized','Position', [0.01 button(4) 0.08 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Specify phase data in NIfTI (.nii or .nii.gz)');
    h.dataIO.edit.inputData1 = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','edit',...
        'units','normalized','position',[0.10 button(4) 0.09 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.dataIO.button.inputData1 = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.2 button(4) 0.05 height],...
        'backgroundcolor','white');
    
    % input data2
    h.dataIO.text.inputData2 = uicontrol(h.StepsPanel.dataIO,...
        'Style','text','String','Magn. data:',...
        'Units','normalized','Position', [0.26 button(4) 0.08 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Specify magnitude data in NIfTI (.nii or .nii.gz)');
    h.dataIO.edit.inputData2 = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','edit',...
        'units','normalized','position',[0.35 button(4) 0.09 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.dataIO.button.inputData2 = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.45 button(4) 0.05 height],...
        'backgroundcolor','white');
    
    % input data3
    h.dataIO.text.inputData3 = uicontrol(h.StepsPanel.dataIO,...
        'Style','text','String','Weights:',...
        'Units','normalized','Position', [0.51 button(4) 0.08 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Specify fieldmapsd data');
    h.dataIO.edit.inputData3 = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','edit',...
        'units','normalized','position',[0.6 button(4) 0.09 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white',...
        'Enable','on');
    h.dataIO.button.inputData3 = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.7 button(4) 0.05 height],...
        'backgroundcolor','white',...
        'Enable','on');
    
    % input header
    h.dataIO.text.inputHeader = uicontrol(h.StepsPanel.dataIO,...
        'Style','text','String','Header:',...
        'Units','normalized','Position', [0.76 button(4) 0.07 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Specify qsm_hub header');
    h.dataIO.edit.inputHeader = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','edit',...
        'units','normalized','position',[0.84 button(4) 0.09 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.dataIO.button.inputHeader = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.94 button(4) 0.05 height],...
        'backgroundcolor','white');

    % output basename
    h.dataIO.text.output = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','text','String','Output basename:',...
        'units','normalized','Position',[0.01 button(3) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Output basename. Default output basename is the input directory followed by /output/squirrel');
    h.dataIO.edit.output = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(3) 0.68 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.dataIO.button.output = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 button(3) 0.1 height],...
        'backgroundcolor','white');

    % invert pahse
    h.dataIO.checkbox.invertPhase = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','checkbox','String','Invert phase data',...
        'units','normalized','Position',[0.01 button(2) 0.3 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Reverse the direction of frequency shift. This could be useful if your previous QSM results show negative values in deep gray matter structure.');

    % brain extraction
    h.dataIO.checkbox.brainExtraction = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','checkbox','String','FSL brain extraction',...
        'units','normalized','Position',[0.33 button(2) 0.3 height],...
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
set(h.dataIO.button.inputData1,         'Callback', {@ButtonOpen_Callback,h,'inputdata1'});
set(h.dataIO.button.inputData2,         'Callback', {@ButtonOpen_Callback,h,'inputdata2'});
set(h.dataIO.button.inputData3,         'Callback', {@ButtonOpen_Callback,h,'inputdata3'});
set(h.dataIO.button.inputHeader,        'Callback', {@ButtonOpen_Callback,h,'header'});
    
end