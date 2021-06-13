%% h = sepia_handle_panel_utility_manage_dependency(hParent,hFig,h,position)
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
% Description: This GUI function creates a panel for the utility function
% 'Get qsm_hub header'
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 June 2021
% Date modified: 
%
%
function h = sepia_handle_panel_utility_manage_dependency(hParent,h,position)

open_icon = imread('folder@0,3x.jpg');
open_icon = imresize(open_icon,[1 1]*16);

%% layout of the panel
nrow        = 10;
rspacing    = 0.03;
ncol        = 1;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Check dependency
FANSI_HOME      = [];
MEDI_HOME       = [];
STISuite_HOME  	= [];
SEGUE_HOME      = [];
ROMEO_HOME      = [];

SpecifyToolboxesDirectory;

%% create panel
% set Parent of all related controls
h.Utility.panel.magageDependency = uipanel(hParent,'Title','Manage Dependencies',...
    'Position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

parent_panel = h.Utility.panel.magageDependency;

wratio = [0.2,0.75,0.05];

    % Dependency 1: FANSI directory input
    [h.Utility.magageDependency.text.FANSIDir,h.Utility.magageDependency.edit.FANSIDir,h.Utility.magageDependency.button.FANSIDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'FANSI directory:',FANSI_HOME,open_icon,[left bottom(1) width height],wratio);
    
    % Dependency 2: MEDI directory input
    [h.Utility.magageDependency.text.MEDIDir,h.Utility.magageDependency.edit.MEDIDir,h.Utility.magageDependency.button.MEDIDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'MEDI directory:',MEDI_HOME,open_icon,[left bottom(2) width height],wratio);
  
    % Dependency 3: STI Suite directory input
    [h.Utility.magageDependency.text.STISuiteDir,h.Utility.magageDependency.edit.STISuiteDir,h.Utility.magageDependency.button.STISuiteDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'STI Suite directory:',STISuite_HOME,open_icon,[left bottom(3) width height],wratio);
    
    % Dependency 4: SEGUE directory input
    [h.Utility.magageDependency.text.SEGUEDir,h.Utility.magageDependency.edit.SEGUEDir,h.Utility.magageDependency.button.SEGUEDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'SEGUE directory:',SEGUE_HOME,open_icon,[left bottom(4) width height],wratio);
    
    % Dependency 5: ROMEO directory input
    [h.Utility.magageDependency.text.ROMEODir,h.Utility.magageDependency.edit.ROMEODir,h.Utility.magageDependency.button.ROMEODir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'ROMEO directory:',ROMEO_HOME,open_icon,[left bottom(5) width height],wratio);
    
    % run
    h.Utility.magageDependency.button.save = uicontrol('Parent',parent_panel,...
        'Style','pushbutton','String','Save',...
        'units','normalized','position',[0.79 bottom(10) 0.2 height],...
        'backgroundcolor','white','enable','off');
    
%% set callback functions
set(h.Utility.magageDependency.button.FANSIDir,  	'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.FANSIDir});
set(h.Utility.magageDependency.button.MEDIDir,      'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.MEDIDir});
set(h.Utility.magageDependency.button.STISuiteDir,	'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.STISuiteDir});
set(h.Utility.magageDependency.button.SEGUEDir,     'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.SEGUEDir});
set(h.Utility.magageDependency.button.ROMEODir,     'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.ROMEODir});
set(h.Utility.magageDependency.button.save,       	'Callback', {@PushbuttonSave_Utility_magageDependency_Callback,h});
end

%% Callback functions
% 'open' button callback
function open_directory_Callback(source,eventdata,h_edit)

% get directory for NIfTI or DICOM files
pathDir = uigetdir;

if pathDir ~= 0
    % set input edit field for display
    set(h_edit, 'String', pathDir);
end

end

function PushbuttonSave_Utility_magageDependency_Callback(source,eventdata,h)

SpecifyToolboxesDirectory;

directory_text = fileread(fullfile(SEPIA_HOME,'SpecifyToolboxesDirectory.m'));


if ~exist('FANSI_HOME','var')   % scenario 1: if such variable doesn't exist yet
    gui_FANSI_HOME = fileparts(get(h.Utility.magageDependency.edit.FANSIDir,'String'));
    % add string to the end of the file
    
elseif isempty(FANSI_HOME)      % scenario 2: if such variable is empty
    gui_FANSI_HOME = fileparts(get(h.Utility.magageDependency.edit.FANSIDir,'String'));
else
    % check if the variable is the same as in the GUI
    curr_FANSI_HOME = fileparts(FANSI_HOME);
    gui_FANSI_HOME  = fileparts(get(h.Utility.magageDependency.edit.FANSIDir,'String'));
    isIdentical     = strcmp(curr_FANSI_HOME,gui_FANSI_HOME);
end




idx  = regexp(directory_text,'FANSI_HOME');

end