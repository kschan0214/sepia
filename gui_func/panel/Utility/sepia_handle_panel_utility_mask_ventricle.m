%% h = sepia_handle_panel_utility_mask_ventricle(hParent,hFig,h,position)
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
% 'Get ventricle mask'
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 25 May 2018
% Date last modified: 12 June 2018
%
%
function h = sepia_handle_panel_utility_mask_ventricle(hParent,h,position)
% set Parent of all related controls
h.Utility.panel.csfMask = uipanel(hParent,'Title','Get ventricle mask',...
    'Position',position,...
    'backgroundcolor',get(h.fig,'color'),...
    'visible','off');
    
    % NIfTI file input
    h.Utility.csfMask.text.niftiInput = uicontrol(h.Utility.panel.csfMask,...
        'Style','text','String','mGRE magnitude',...
        'Units','normalized','Position', [0.01 0.85 0.2 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Multi-echo GRE (4D) magnitude NIfTI file.');
    h.Utility.csfMask.edit.niftiInput = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.85 0.68 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.csfMask.button.niftiInput = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 0.85 0.1 0.14],...
        'backgroundcolor','white');
    
    % mask file input
    h.Utility.csfMask.text.maskInput = uicontrol(h.Utility.panel.csfMask,...
        'Style','text','String','Brain mask',...
        'Units','normalized','Position', [0.01 0.70 0.2 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Brain mask NIfTI file (optional)');
    h.Utility.csfMask.edit.maskInput = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.70 0.68 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.csfMask.button.maskInput = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 0.70 0.1 0.14],...
        'backgroundcolor','white');
    
    % TE
    % TE file
    h.Utility.csfMask.text.teFile = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','text','String','TE file',...
        'units','normalized','Position',[0.01 0.55 0.2 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Can be either qsm_hub header (.mat) or MRIConvert text (.txt) file');
    h.Utility.csfMask.edit.teFile = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.55 0.68 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.csfMask.button.teFile = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 0.55 0.1 0.14],...
        'backgroundcolor','white');
    % user TE input
    h.Utility.csfMask.text.userTE = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','text','String','or user input TEs (s)',...
        'units','normalized','Position',[0.01 0.4 0.2 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Specify the echo time (TE, in s) of each echo');
    h.Utility.csfMask.edit.userTE = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.4 0.78 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white',...
        'String','[ ]');
    
    % output directory
    h.Utility.csfMask.text.outputDir = uicontrol(h.Utility.panel.csfMask,...
        'Style','text','String','Output dir:',...
        'Units','normalized','Position', [0.01 0.25 0.2 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Output directory where the output mask file will be stored');
    h.Utility.csfMask.edit.outputDir = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','edit',...
        'units','normalized','position',[0.21 0.25 0.68 0.14],...
        'HorizontalAlignment','left',...
        'backgroundcolor','white');
    h.Utility.csfMask.button.outputDir = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','pushbutton','String','open',...
        'units','normalized','position',[0.89 0.25 0.1 0.14],...
        'backgroundcolor','white');
    
    % run
    h.Utility.csfMask.button.run = uicontrol('Parent',h.Utility.panel.csfMask,...
        'Style','pushbutton','String','Run',...
        'units','normalized','position',[0.89 0.05 0.1 0.18],...
        'backgroundcolor','white');

%% set callbacks
set(h.Utility.csfMask.button.teFile,        'Callback',             {@ButtonOpen_Utility_csfMask_Callback,h,'te'});
set(h.Utility.csfMask.button.niftiInput,    'Callback',             {@ButtonOpen_Utility_csfMask_Callback,h,'nitfi'});
set(h.Utility.csfMask.button.maskInput,     'Callback',             {@ButtonOpen_Utility_csfMask_Callback,h,'mask'});
set(h.Utility.csfMask.button.outputDir,     'Callback',             {@ButtonOpen_Utility_csfMask_Callback,h,'output'});
set(h.Utility.csfMask.button.run,           'Callback',             {@PushbuttonRun_Utility_csfMask_Callback,h});
    
end

%% Callback functions
function ButtonOpen_Utility_csfMask_Callback(source,eventdata,h,field)
% get input file/directory for getHeader utility function
% global h

switch field
    case 'te'
        % te file can be text file or mat file
        [tefileName,pathDir] = uigetfile({'*.mat;*.txt'},'Select TE file');
        
        % display file directory
        if pathDir ~= 0
            set(h.Utility.csfMask.edit.teFile,'String',fullfile(pathDir,tefileName));
        end

    case 'nitfi'
        % read NIfTI file 
        [nitfiName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select NIfTI file');

        if pathDir ~= 0
            set(h.Utility.csfMask.edit.niftiInput,    'String',fullfile(pathDir,nitfiName));
            % automatically set default output field
            set(h.Utility.csfMask.edit.outputDir,     'String',[pathDir filesep 'output']);
        end
        
    case 'mask'
        % read NIfTI file 
        [nitfiName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select mask file');

        if pathDir ~= 0
            set(h.Utility.csfMask.edit.maskInput,    'String',fullfile(pathDir,nitfiName));
        end
        
    case 'output'
        
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.Utility.csfMask.edit.outputDir,     'String',pathDir);
        end
end

end

function PushbuttonRun_Utility_csfMask_Callback(source,event,h)
% Callback function to get lateral ventricle mask

% global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

try

% get NIfTI file (if any)
niftiFile   = get(h.Utility.csfMask.edit.niftiInput, 'String');
% get output directory, assume the same directory as input directory/file
outputDir = get(h.Utility.csfMask.edit.outputDir, 'String');
% mask file
maskFullName    = get(h.Utility.csfMask.edit.maskInput,'String');

if isempty(niftiFile)
    % This function requires input file/directory
    error('Please specify a NIfTI file.');
else

    % try to get TE variable in this stage
    teFullName = get(h.Utility.csfMask.edit.teFile, 'String');
    teUserInput = get(h.Utility.csfMask.edit.userTE, 'String');
    if ~isempty(teFullName)
        % get data type
        [~,~,ext] = fileparts(teFullName);

        if strcmpi(ext,'.mat')
            % if mat file then try to load 'TE' directly
            try load(teFullName,'TE');  catch; error('No variable named TE.'); end
        else
            % if text file the try to read the TEs line by line
            TE = readTEfromText(teFullName);
            TE = TE(:);
        end
    else
        % read user input array
        TE = str2num(teUserInput);
        TE = TE(:);
    end
    % in case TE cannot be found
    if isempty(TE) || isnan(TE(1))
        error('Incorrect TE format.');
    end
    
    GetVentricleMaskIOWrapper(niftiFile,outputDir,maskFullName,TE);
    
end

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end
    

% Disable the pushbutton to prevent doubel click
set(source,'Enable','on');

end