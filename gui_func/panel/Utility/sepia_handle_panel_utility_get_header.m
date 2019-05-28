%% h = sepia_handle_panel_utility_get_header(hParent,hFig,h,position)
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
% Date created: 24 May 2018
% Date modified: 12 June 2018
% Date modified: 24 May 2019
%
%
function h = sepia_handle_panel_utility_get_header(hParent,h,position)

% define maximum level of options and spacing between options
nlevel = 10;
spacing = 0.02;
height = (1-(nlevel+1)*spacing)/nlevel;
button = (height+spacing:height+spacing:(height+spacing)*nlevel) - height;


% set Parent of all related controls
h.Utility.panel.getHeader = uipanel(hParent,'Title','Get Sepia header',...
    'Position',position,...
    'backgroundcolor',get(h.fig,'color'));

    % Option 1: DICOM directory input
    h.Utility.getHeader.text.dicomInput = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','Op 1: Input DICOM dir:',...
        'Units','normalized','Position', [0.01 button(10) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',...
        strcat('Input directory contains all (magnitude and phase) mGRE DICOM files. ',...
         'No user input will be needed if your data input is a DICOM directory.'));
    h.Utility.getHeader.edit.dicomInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(10) 0.68 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.dicomInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 button(10) 0.1 height],...
        'backgroundcolor','white');
    
    % Option 2: NIfTI directory input
    h.Utility.getHeader.text.niftiDirInput = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','Op 2: Input NIfTI dir:',...
        'Units','normalized','Position', [0.01 button(9) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',...
        strcat('Input directory contains ONLY NIfTI and JSON/text files. ',...
         'No user input will be needed if your data input is a directory.'));
    h.Utility.getHeader.edit.niftiDirInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(9) 0.68 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.niftiDirInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 button(9) 0.1 height],...
        'backgroundcolor','white');
    
    % Option 3: NIfTI file input
    h.Utility.getHeader.text.niftiInput = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','Op 3: Select an NIfTI file:',...
        'Units','normalized','Position', [0.01 button(8) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Input NIfTI file of mGRE data. User input is allowed');
    h.Utility.getHeader.edit.niftiInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(8) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.niftiInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.41 button(8) 0.08 height],...
        'backgroundcolor','white');
    % TE file
    h.Utility.getHeader.text.teFile = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','text','String','and TE file(s)',...
        'units','normalized','Position',[0.51 button(8) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Can be either Sepia header (.mat), MRIConvert text (.txt) or dicm2nii or dcm2niix JSON file(s)');
    h.Utility.getHeader.edit.teFile = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.71 button(8) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.teFile = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.91 button(8) 0.08 height],...
        'backgroundcolor','white');
    
    % output directory
    h.Utility.getHeader.text.outputDir = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','Output basename:',...
        'Units','normalized','Position', [0.01 button(7) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Output directory where the MAT header file will be stored');
    h.Utility.getHeader.edit.outputDir = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(7) 0.68 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.outputDir = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 button(7) 0.1 height],...
        'backgroundcolor','white');
    
    
    % User input
    h.Utility.getHeader.text.user = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','User defined input. These values will overide the information detected from the input data.',...
        'Units','normalized','Position', [0.01 button(5) 0.9 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',...
        strcat('User input will replace the header information extracted from the input NIfTI header. ',...
               'Information about echo times has to be provided along with the NIfTI file.'));
    % B0 strength
    h.Utility.getHeader.text.userB0 = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','B0 strength (T)',...
        'Units','normalized','Position', [0.01 button(4) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Magnetic field strength in Tesla (default: 3).');
    h.Utility.getHeader.edit.userB0 = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'String','',...
        'units','normalized','position',[0.21 button(4) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    % B0 direction
    h.Utility.getHeader.text.userB0dir = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','B0 direction [x,y,z]',...
        'Units','normalized','Position', [0.49 button(4) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Magnetic field direction');
    h.Utility.getHeader.edit.userB0dir = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit','String','[ ]',...
        'units','normalized','position',[0.71 button(4) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    
    % voxel size
    h.Utility.getHeader.text.userVoxelSize = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','Voxel size (x,y,z) (mm)',...
        'Units','normalized','Position', [0.01 button(3) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Voxel size of image in order of [x,y,z]');
    h.Utility.getHeader.edit.userVoxelSize = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit','String','[ ]',...
        'units','normalized','position',[0.21 button(3) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    
    % TE
    % user TE input
    h.Utility.getHeader.text.userTE = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','text','String','or user input TEs (s)',...
        'units','normalized','Position',[0.01 button(2) 0.2 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Specify the echo time (TE, in s) of each echo');
    h.Utility.getHeader.edit.userTE = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 button(2) 0.78 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white',...
        'String','[ ]');
    
    
    
    % run
    h.Utility.getHeader.button.run = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','Save header',...
        'units','normalized','position',[0.79 button(1) 0.2 height],...
        'backgroundcolor','white');
    
%% set callback functions
set(h.Utility.getHeader.button.teFile,          'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'te'});
set(h.Utility.getHeader.button.dicomInput,      'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'dicom'});
set(h.Utility.getHeader.button.niftiInput,      'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'nitfi'});
set(h.Utility.getHeader.button.niftiDirInput,   'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'niftidir'});
set(h.Utility.getHeader.button.outputDir,       'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'output'});
set(h.Utility.getHeader.button.run,             'Callback', {@PushbuttonRun_Utility_getHeader_Callback,h});
end

%% Callback functions
% 'open' button callback
function ButtonOpen_Utility_getHeader_Callback(source,eventdata,h,field)
% get input file/directory for getHeader utility function
% global h

prefix = 'sepia';

switch field
    case 'te'
        % te file can be text file or mat file
        [tefileName,pathDir] = uigetfile({'*.txt;*.mat;*.json'},'Select TE file(s)','MultiSelect', 'on');
        
        fullname = [];
        if iscell(tefileName)
            for kf = 1:length(tefileName)
                fullname = [fullname, fullfile(pathDir,tefileName{kf}) ';'];
            end
        else
            fullname = fullfile(pathDir,tefileName);
        end
        
        % display file directory
        if pathDir ~= 0
            set(h.Utility.getHeader.edit.teFile,'String',fullname);
        end

    case 'nitfi'
        % read NIfTI file 
        [nitfiName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select mask file');
        
        if pathDir ~= 0
            
            input_nii = load_untouch_nii(fullfile(pathDir,nitfiName));
            [~,B0_dir,voxelSize,~,~,~,~] = SyntheticQSMHubHeader(input_nii);
            b0_dir_str = sprintf('[%.4f, %.4f, %.4f]',B0_dir(1),B0_dir(2),B0_dir(3));
            voxel_size_str = sprintf('[%.4f, %.4f, %.4f]',voxelSize(1),voxelSize(2),voxelSize(3));
            set(h.Utility.getHeader.edit.userB0dir,     'String',b0_dir_str);
            set(h.Utility.getHeader.edit.userVoxelSize, 'String',voxel_size_str);
            
            set(h.Utility.getHeader.edit.niftiInput,    'String',fullfile(pathDir,nitfiName));
            % automatically set default output field
            set(h.Utility.getHeader.edit.outputDir,     'String',[pathDir prefix]);
            % empty DICOM input field
            set(h.Utility.getHeader.edit.dicomInput,    'String',[]);
            % empty NIfTi directory field
            set(h.Utility.getHeader.edit.niftiDirInput, 'String',[]);
        end
        
    case 'dicom'
        % get DICOM directory
        pathDir = uigetdir;

        if pathDir ~= 0
            % set input edit field for display
            set(h.Utility.getHeader.edit.dicomInput,    'String',pathDir);
            % automatically set default output field
            set(h.Utility.getHeader.edit.outputDir,     'String',fullfile(pathDir,prefix));
            % empty NIfTI input field
            set(h.Utility.getHeader.edit.niftiInput,    'String',[]);
            % empty NIfTi directory field
            set(h.Utility.getHeader.edit.niftiDirInput, 'String',[]);
            % empty TE files field
            set(h.Utility.getHeader.edit.teFile,        'String',[]);
            % empty user defined field
            set(h.Utility.getHeader.edit.userB0,        'String',[]);
            set(h.Utility.getHeader.edit.userB0dir,     'String','[ ]');
            set(h.Utility.getHeader.edit.userVoxelSize, 'String','[ ]');
            set(h.Utility.getHeader.edit.userTE,        'String','[ ]');
        end
        
    case 'niftidir'
        % get directory
        pathDir = uigetdir;

        if pathDir ~= 0
            % set input edit field for display
            set(h.Utility.getHeader.edit.niftiDirInput, 'String',pathDir);
            % automatically set default output field
            set(h.Utility.getHeader.edit.outputDir,     'String',fullfile(pathDir,prefix));
            % empty NIfTI input field
            set(h.Utility.getHeader.edit.niftiInput,    'String',[]);
            % empty DICOM directory field
            set(h.Utility.getHeader.edit.dicomInput,    'String',[]);
            % empty TE files field
            set(h.Utility.getHeader.edit.teFile,        'String',[]);
            % empty user defined field
            set(h.Utility.getHeader.edit.userB0,        'String',[]);
            set(h.Utility.getHeader.edit.userB0dir,     'String','[ ]');
            set(h.Utility.getHeader.edit.userVoxelSize, 'String','[ ]');
            set(h.Utility.getHeader.edit.userTE,        'String','[ ]');

        end
        
    case 'output'
        
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.Utility.getHeader.edit.outputDir,     'String',[pathDir filesep prefix]);
        end
end

end

function PushbuttonRun_Utility_getHeader_Callback(source,eventdata,h)
% Callback function to detect and save header functino for sepia

% global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

try 
sepia_addpath('dicom');

% input files/directory
% get DICOM directory (if any)
dicomDir    = get(h.Utility.getHeader.edit.dicomInput,'String'); 
% get Nifti directory (if any)
niftiDir    = get(h.Utility.getHeader.edit.niftiDirInput,'String'); 
% get output directory, assume the same directory as input directory/file
outputDir = get(h.Utility.getHeader.edit.outputDir, 'String');
% get NIfTI file (if any)
niftiFile   = get(h.Utility.getHeader.edit.niftiInput, 'String');
% get TE files (if any)
teFullName = get(h.Utility.getHeader.edit.teFile, 'String');
% re-organise TE file name variable to cell
if ~isempty(teFullName)
    kInd = strfind(teFullName,';');
    if kInd ~= 0
        for kf = 1:length(kInd)
            if kf==1
                tmpName{kf} = teFullName(1:kInd(kf)-1);
            else
                tmpName{kf} = teFullName(kInd(kf-1)+1:kInd(kf)-1);
            end
        end
        teFullName = tmpName;
    else
        tmpName{1} = teFullName;
        teFullName = tmpName;
    end
end

% user defined input
userDefine = [];
% get user input magnetic field strength
B0           = str2double(get(h.Utility.getHeader.edit.userB0, 'String'));
% str2double returns NaN if input field is empty, replace it by []
if ~isnan(B0)
    userDefine.B0 = B0;
end
% get user input magnetic field direction
B0_dir       = str2num(get(h.Utility.getHeader.edit.userB0dir, 'String'));
% check validity of input B0 direction
if length(B0_dir) ~= 3 && ~isempty(B0_dir) && isempty(dicomDir)
    error('The B0 direction has 3 elements [x,y,z]');
end
if ~isempty(B0_dir)
    userDefine.B0_dir = B0_dir;
end
% get user input voxel size
voxelSize    = str2num(get(h.Utility.getHeader.edit.userVoxelSize, 'String'));
% check validity of input voxel size
if length(voxelSize) ~= 3 && ~isempty(voxelSize) && isempty(dicomDir)
    error(['Sepia currently works with 3D data only. ' 
           'Please specify the voxel size in all three dimensions']);
end
if ~isempty(voxelSize)
    userDefine.voxelSize = voxelSize;
end
% get user input TE
TE           = str2num(get(h.Utility.getHeader.edit.userTE, 'String'));
if ~isempty(TE)
    userDefine.TE = TE;
end

% if no dicom directory detected then get te from user input
if ~isempty(dicomDir)
    
    save_sepia_header(dicomDir,userDefine,outputDir);
    
elseif ~isempty(niftiDir)
    
    save_sepia_header(niftiDir,userDefine,outputDir);
    
elseif ~isempty(niftiFile)
    
    input.nifti = niftiFile;
    input.TEFileList = teFullName;
    
    save_sepia_header(input,userDefine,outputDir);
    
else

    % This function requires input file/directory
    error('Please specify a directory containing DICOM/NIfTI files or a NIfTI file.');
    
end

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end

disp('SEPIA header is saved!');

% re-enable the pushbutton 
set(source,'Enable','on');

end