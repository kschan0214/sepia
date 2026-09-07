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
% Date created: 22 June 2021
% Date modified: 
%
%
function h = sepia_handle_panel_utility_manage_dependency(hParent,h,position)

open_icon = sepia_theme_open_icon();

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
MRITOOLS_HOME   = [];
MRISC_HOME      = [];
ANTS_HOME       = [];
HEIDI_HOME      = [];
ChiSepNet_HOME  = [];

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
        'FANSI Home:',FANSI_HOME,open_icon,[left bottom(1) width height],wratio);
    
    % Dependency 2: MEDI directory input
    [h.Utility.magageDependency.text.MEDIDir,h.Utility.magageDependency.edit.MEDIDir,h.Utility.magageDependency.button.MEDIDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'MEDI Home:',MEDI_HOME,open_icon,[left bottom(2) width height],wratio);
  
    % Dependency 3: STI Suite directory input
    [h.Utility.magageDependency.text.STISuiteDir,h.Utility.magageDependency.edit.STISuiteDir,h.Utility.magageDependency.button.STISuiteDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'STI Suite Home:',STISuite_HOME,open_icon,[left bottom(3) width height],wratio);
    
    % Dependency 4: SEGUE directory input
    [h.Utility.magageDependency.text.SEGUEDir,h.Utility.magageDependency.edit.SEGUEDir,h.Utility.magageDependency.button.SEGUEDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'SEGUE Home:',SEGUE_HOME,open_icon,[left bottom(4) width height],wratio);
    
    % Dependency 5: ROMEO directory input
    [h.Utility.magageDependency.text.MRITOOLSDir,h.Utility.magageDependency.edit.MRITOOLSDir,h.Utility.magageDependency.button.MRITOOLSDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'MRITOOLS Home:',MRITOOLS_HOME,open_icon,[left bottom(5) width height],wratio);

    % Dependency 6: ROMEO directory input
    [h.Utility.magageDependency.text.MRISuscCalcDir,h.Utility.magageDependency.edit.MRISuscCalcDir,h.Utility.magageDependency.button.MRISuscCalcDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'MRI susc. calc. Home:',MRISC_HOME,open_icon,[left bottom(6) width height],wratio);
    
    % Dependency ANTs
    [h.Utility.magageDependency.text.ANTsDir,h.Utility.magageDependency.edit.ANTsDir,h.Utility.magageDependency.button.ANTsDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'ANTs Home:',ANTS_HOME,open_icon,[left bottom(7) width height],wratio);

    % Dependency HEIDI
    [h.Utility.magageDependency.text.HEIDIDir,h.Utility.magageDependency.edit.HEIDIDir,h.Utility.magageDependency.button.HEIDIDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'HEIDI Home:',HEIDI_HOME,open_icon,[left bottom(8) width height],wratio);

    % Dependency Chi-separation
    [h.Utility.magageDependency.text.ChiSepNetDir,h.Utility.magageDependency.edit.ChiSepNetDir,h.Utility.magageDependency.button.ChiSepNetDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'Chi-separation Home:',ChiSepNet_HOME,open_icon,[left bottom(9) width height],wratio);

    % run
    h.Utility.magageDependency.button.save = uicontrol('Parent',parent_panel,...
        'Style','pushbutton','String','Save',...
        'units','normalized','position',[0.79 bottom(10) 0.2 height],...
        'backgroundcolor',sepia_theme_bg_color(),'foregroundcolor',sepia_theme_fg_color(),'enable','on');

%% set callback functions
set(h.Utility.magageDependency.button.FANSIDir,         'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.FANSIDir});
set(h.Utility.magageDependency.button.MEDIDir,          'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.MEDIDir});
set(h.Utility.magageDependency.button.STISuiteDir,      'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.STISuiteDir});
set(h.Utility.magageDependency.button.SEGUEDir,         'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.SEGUEDir});
set(h.Utility.magageDependency.button.MRITOOLSDir,      'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.MRITOOLSDir});
set(h.Utility.magageDependency.button.MRISuscCalcDir,   'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.MRISuscCalcDir});
set(h.Utility.magageDependency.button.ANTsDir,       	'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.ANTsDir});
set(h.Utility.magageDependency.button.HEIDIDir,       	'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.HEIDIDir});
set(h.Utility.magageDependency.button.ChiSepNetDir,   	'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.ChiSepNetDir});
set(h.Utility.magageDependency.button.save,             'Callback', {@PushbuttonSave_Utility_magageDependency_Callback,h});
end

%% Callback functions
% 'open' button callback
function open_directory_Callback(source,eventdata,h_edit)

% get directory for NIfTI or DICOM files
pathDir = uigetdir;

if pathDir ~= 0
    % set input edit field for display
    set(h_edit, 'String', [pathDir filesep]);
end

end

function PushbuttonSave_Utility_magageDependency_Callback(source,eventdata,h)

dependency_homes = {'FANSI_HOME','MEDI_HOME','STISuite_HOME','SEGUE_HOME','MRITOOLS_HOME','MRISC_HOME', 'ANTS_HOME', 'HEIDI_HOME', 'ChiSepNet_HOME'};
gui_handles      = {'FANSIDir'  ,'MEDIDir'  ,'STISuiteDir'  ,'SEGUEDir'  ,'MRITOOLSDir'  ,'MRISuscCalcDir', 'ANTsDir', 'HEIDIDir', 'ChiSepNetDir'};

sepia_universal_variables;

isOverWrite = false;

for k = 1:length(gui_handles)

% get string from GUI
gui_field = get(h.Utility.magageDependency.edit.(gui_handles{k}),'String');

% if GUI is not empty, then allows changes
if ~isempty( gui_field )
    isUpdated   = update_toolbox_directory_entry(SEPIA_HOME, dependency_homes{k}, gui_field);
    isOverWrite = isOverWrite || isUpdated;
end
end

if isOverWrite
    disp('The paths are save in SpecifyToolboxesDirectory.m!')
else
    disp('No changes to SpecifyToolboxesDirectory.m were needed.')
end

end