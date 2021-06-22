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
MRISC_HOME      = [];

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
    [h.Utility.magageDependency.text.ROMEODir,h.Utility.magageDependency.edit.ROMEODir,h.Utility.magageDependency.button.ROMEODir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'ROMEO Home:',ROMEO_HOME,open_icon,[left bottom(5) width height],wratio);

    % Dependency 6: ROMEO directory input
    [h.Utility.magageDependency.text.MRISuscCalcDir,h.Utility.magageDependency.edit.MRISuscCalcDir,h.Utility.magageDependency.button.MRISuscCalcDir] = ...
        sepia_construct_text_edit_button(parent_panel,...
        'MRI susc. calc. Home:',MRISC_HOME,open_icon,[left bottom(6) width height],wratio);
    
    % run
    h.Utility.magageDependency.button.save = uicontrol('Parent',parent_panel,...
        'Style','pushbutton','String','Save',...
        'units','normalized','position',[0.79 bottom(10) 0.2 height],...
        'backgroundcolor','white','enable','on');
    
%% set callback functions
set(h.Utility.magageDependency.button.FANSIDir,         'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.FANSIDir});
set(h.Utility.magageDependency.button.MEDIDir,          'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.MEDIDir});
set(h.Utility.magageDependency.button.STISuiteDir,      'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.STISuiteDir});
set(h.Utility.magageDependency.button.SEGUEDir,         'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.SEGUEDir});
set(h.Utility.magageDependency.button.ROMEODir,         'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.ROMEODir});
set(h.Utility.magageDependency.button.MRISuscCalcDir,   'Callback', {@open_directory_Callback,h.Utility.magageDependency.edit.MRISuscCalcDir});
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

dependency_homes = {'FANSI_HOME','MEDI_HOME','STISuite_HOME','SEGUE_HOME','ROMEO_HOME','MRISC_HOME'};
gui_handles      = {'FANSIDir'  ,'MEDIDir'  ,'STISuiteDir'  ,'SEGUEDir'  ,'ROMEODir'  ,'MRISuscCalcDir'};

SpecifyToolboxesDirectory;

% get all the text from SpecifyToolboxesDirectory.m 
fid             = fopen( fullfile(SEPIA_HOME,'SpecifyToolboxesDirectory.m') );
directory_text  = textscan( fid, '%s', 'Delimiter','\n', 'CollectOutput',true );
fclose( fid );

isOverWrite = false;

for k = 1:length(gui_handles)

% get string from GUI
gui_field = get(h.Utility.magageDependency.edit.(gui_handles{k}),'String');

% if GUI is not empty, then allows changes
if ~isempty( gui_field )
    
    gui_HOME = fileparts(fullfile(gui_field,filesep));
    
    % default update is false
    isUpdateHome = false;
    
    % check if changing is needed for the following conditions
    if ~exist(dependency_homes{k},'var')                % scenario 1: if such variable doesn't exist yet
        isUpdateHome = true;

    elseif isempty(eval(dependency_homes{k}))           % scenario 2: if such variable is empty
        isUpdateHome = true;
    else                                                % scenario 3: check if the variable is the same as in the GUI
        curr_HOME = fileparts(eval(dependency_homes{k}));
        isUpdateHome = ~strcmp(curr_HOME,gui_HOME);    % if not identical then update
    end
    
    % update SpecifyToolboxesDirectory.m
    if isUpdateHome
        
        % check if the file contains the variable name that is about to be changed 
        % if so, and if the 1st char is not '%' then comment the line out
        for j = 1:length(directory_text{1})
            
            isContain = ContainName(directory_text{1}{j},lower(dependency_homes{k}));
            
            if isContain && ~strcmp(directory_text{1}{j}(1),'%')
                % insert a '%' to comment the line out
                directory_text{1}{j} = ['% ' directory_text{1}{j}];
            end
                
        end
        % insert the variable to the end of the file
        directory_text{1}{j+1} = sprintf('%s = ''%s'';',dependency_homes{k}, fullfile(gui_HOME,filesep));
        
        isOverWrite = true;
    end
    
end
end

% overwrite SpecifyToolboxesDirectory.m
if isOverWrite
    
    fid = fopen( fullfile(SEPIA_HOME,'SpecifyToolboxesDirectory.m'), 'w');
    for j = 1:length(directory_text{1})
        fprintf( fid, '%s\n', directory_text{1}{j} );
    end
    fclose( fid );

end


end