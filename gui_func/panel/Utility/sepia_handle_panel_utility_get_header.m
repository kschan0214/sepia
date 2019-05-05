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
% Date last modified: 12 June 2018
%
%
function h = sepia_handle_panel_utility_get_header(hParent,h,position)
% set Parent of all related controls
h.Utility.panel.getHeader = uipanel(hParent,'Title','Get Sepia header',...
    'Position',position,...
    'backgroundcolor',get(h.fig,'color'));

    % DICOM directory input
    h.Utility.getHeader.text.dicomInput = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','Input DICOM dir:',...
        'Units','normalized','Position', [0.01 0.88 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',...
        strcat('Input directory contains all (magnitude and phase) mGRE DICOM files. ',...
         'No user input will be needed if your data input is a DICOM directory.'));
    h.Utility.getHeader.edit.dicomInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.88 0.68 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.dicomInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 0.88 0.1 0.1],...
        'backgroundcolor','white');
    
    % NIfTI file input
    h.Utility.getHeader.text.niftiInput = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','or input NIfTI file:',...
        'Units','normalized','Position', [0.01 0.77 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Input NIfTI file of mGRE data. User input is allowed');
    h.Utility.getHeader.edit.niftiInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.77 0.68 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.niftiInput = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 0.77 0.1 0.1],...
        'backgroundcolor','white');
    
    % User input
    h.Utility.getHeader.text.user = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','or/and user-defined parameter(s)',...
        'Units','normalized','Position', [0.01 0.66 0.9 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',...
        strcat('User input will replace the header information extracted from the input NIfTI header. ',...
               'Information about echo times has to be provided along with the NIfTI file.'));
    % B0 strength
    h.Utility.getHeader.text.userB0 = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','B0 strength (T)',...
        'Units','normalized','Position', [0.01 0.55 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Magnetic field strength in Tesla (default: 3).');
    h.Utility.getHeader.edit.userB0 = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'String','3',...
        'units','normalized','position',[0.21 0.55 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    % B0 direction
    h.Utility.getHeader.text.userB0dir = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','B0 direction [x,y,z]',...
        'Units','normalized','Position', [0.49 0.55 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Magnetic field direction');
    h.Utility.getHeader.edit.userB0dir = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit','String','[ ]',...
        'units','normalized','position',[0.71 0.55 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    
    % voxel size
    h.Utility.getHeader.text.userVoxelSize = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','Voxel size (x,y,z) (mm)',...
        'Units','normalized','Position', [0.01 0.44 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Voxel size of image in order of [x,y,z]');
    h.Utility.getHeader.edit.userVoxelSize = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit','String','[ ]',...
        'units','normalized','position',[0.21 0.44 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    
    % TE
    % TE file
    h.Utility.getHeader.text.teFile = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','text','String','TE file',...
        'units','normalized','Position',[0.01 0.33 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Can be either Sepia header (.mat), MRIConvert text (.txt) or dicm2nii or dcm2niix JSON file(s)');
    h.Utility.getHeader.edit.teFile = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.33 0.68 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.teFile = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 0.33 0.1 0.1],...
        'backgroundcolor','white');
    % user TE input
    h.Utility.getHeader.text.userTE = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','text','String','or user input TEs (s)',...
        'units','normalized','Position',[0.01 0.22 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Specify the echo time (TE, in s) of each echo');
    h.Utility.getHeader.edit.userTE = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.22 0.78 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white',...
        'String','[ ]');
    
    % output directory
    h.Utility.getHeader.text.outputDir = uicontrol(h.Utility.panel.getHeader,...
        'Style','text','String','Output basename:',...
        'Units','normalized','Position', [0.01 0.11 0.2 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Output directory where the MAT header file will be stored');
    h.Utility.getHeader.edit.outputDir = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.11 0.68 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.getHeader.button.outputDir = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 0.11 0.1 0.1],...
        'backgroundcolor','white');
    
    % run
    h.Utility.getHeader.button.run = uicontrol('Parent',h.Utility.panel.getHeader,...
        'Style','pushbutton','String','Save header',...
        'units','normalized','position',[0.79 0.005 0.2 0.1],...
        'backgroundcolor','white');
    
%% set callback functions
set(h.Utility.getHeader.button.teFile,      'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'te'});
set(h.Utility.getHeader.button.dicomInput,  'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'dicom'});
set(h.Utility.getHeader.button.niftiInput,  'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'nitfi'});
set(h.Utility.getHeader.button.outputDir,   'Callback', {@ButtonOpen_Utility_getHeader_Callback,h,'output'});
set(h.Utility.getHeader.button.run,         'Callback', {@PushbuttonRun_Utility_getHeader_Callback,h});
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
        end
        
    case 'dicom'
        % get DICOM directory
        pathDir = uigetdir;

        if pathDir ~= 0
            % set input edit field for display
            set(h.Utility.getHeader.edit.dicomInput,    'String',pathDir);
            % automatically set default output field
            set(h.Utility.getHeader.edit.outputDir,     'String',pathDir);
            % empty NIfTI input field
            set(h.Utility.getHeader.edit.niftiInput,    'String',[]);
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
% Callback function to detect and save header functino for qsm_hub

% global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

try 
sepia_addpath('dicom');
% get DICOM directory (if any)
dicomDir    = get(h.Utility.getHeader.edit.dicomInput,'String'); 
% get output directory, assume the same directory as input directory/file
outputDir = get(h.Utility.getHeader.edit.outputDir, 'String');


% if no dicom directory detected then get te from user input
if isempty(dicomDir)
    
    % get NIfTI file (if any)
    niftiFile   = get(h.Utility.getHeader.edit.niftiInput, 'String');
    
    if isempty(niftiFile)
        % This function requires input file/directory
        error('Please specify a DICOM directory or a NIfTI file.');
    else
        
        % get user input magnetic field strength
        b0          = str2double(get(h.Utility.getHeader.edit.userB0, 'String'));
        % get user input magnetic field direction
        b0dir       = str2num(get(h.Utility.getHeader.edit.userB0dir, 'String'));
        % get user input voxel size
        voxelSize   = str2num(get(h.Utility.getHeader.edit.userVoxelSize, 'String'));

        % check validity of input voxel size
        if length(voxelSize) ~= 3 && ~isempty(voxelSize) && isempty(dicomDir)
            error(['Sepia currently works with 3D data only. ' 
                   'Please specify the voxel size in all three dimensions']);
        end

        % check validity of input B0 direction
        if length(b0dir) ~= 3 && ~isempty(b0dir) && isempty(dicomDir)
            error('The B0 direction has 3 elements [x,y,z]');
        end

        % str2double returns NaN if input field is empty, replace it by []
        if isnan(b0)
            b0=[];
        end

        % try to get TE variable in this stage
        teFullName = get(h.Utility.getHeader.edit.teFile, 'String');
        teUserInput = get(h.Utility.getHeader.edit.userTE, 'String');
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
                [~,~,ext] = fileparts(teFullName{1});
            else
                % get data type
                [~,~,ext] = fileparts(teFullName);
            end

            if strcmpi(ext,'.mat')
                % if mat file then try to load 'TE' directly
                try load(teFullName,'TE');  catch; error('No variable named TE.'); end
            elseif strcmpi(ext, '.txt')
                % if text file the try to read the TEs line by line
                TE = readTE_MRIConvert_Text(teFullName);
            elseif strcmpi(ext, '.json')
                % JSON file(s)
                TE = readTE_JSON(teFullName);
            end
            TE = TE(:).';
        else
            % read user input array
            TE = str2num(teUserInput);
            TE = TE(:);
        end
        % in case TE cannot be found
        if isempty(TE) || isnan(TE(1))
            error('Incorrect TE format.');
        end
        
        % specify the input is NIfTI file
        ExportQMSHubHeaderIOWrapper(niftiFile,outputDir,...
            b0,b0dir,voxelSize,TE,'nifti');
    end
else
    
    % specify the input is DICOM directory
    ExportQMSHubHeaderIOWrapper(dicomDir,outputDir,[],[],[],[],'dicom');
    
end

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end

disp('Done!');

% re-enable the pushbutton 
set(source,'Enable','on');

end