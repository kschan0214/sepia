function ButtonOpen_Callback(source,eventdata,h,field)
% get directory and display it on GUI

% global h

% default output base name
prefix = 'Sepia';

switch field
    case 'mask'
        % only read NIfTI file for mask
        [maskfileName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select mask file');

        if pathDir ~= 0
            set(h.dataIO.edit.maskdir,'String',fullfile(pathDir,maskfileName));
        end
        
    case 'input'
        % get directory for NIfTI or DICOM files
        pathDir = uigetdir;

        if pathDir ~= 0
            % set input edit field for display
            set(h.dataIO.edit.input,    'String',pathDir);
            % automatically set default output field
            set(h.dataIO.edit.output,   'String',[pathDir filesep 'output' filesep prefix]);
            
            % set input data fields empty
            set(h.dataIO.edit.inputData1,    'String',[]);
            set(h.dataIO.edit.inputData2,    'String',[]);
            set(h.dataIO.edit.inputData3,    'String',[]);
            set(h.dataIO.edit.inputHeader,   'String',[]);
        end
        
    case 'output'
        
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.dataIO.edit.output,'String',[pathDir filesep prefix]);
        end
        
    case 'inputdata1'
        % only read NIfTI file for mask
        [fileName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select a NIfTI file');

        if pathDir ~= 0
            % set input edit field for display
            set(h.dataIO.edit.inputData1,    'String',fullfile(pathDir,fileName));
            % automatically set default output field
            set(h.dataIO.edit.output,   'String',[pathDir 'output' filesep prefix]);
            % set input directory field empty
            set(h.dataIO.edit.input,    'String',[]);
        end
        
    case 'inputdata2'
        % only read NIfTI file for mask
        [fileName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select a NIfTI file');

        if pathDir ~= 0
            % set input edit field for display
            set(h.dataIO.edit.inputData2,    'String',fullfile(pathDir,fileName));
            
            % set input directory field empty
            set(h.dataIO.edit.input,    'String',[]);

        end
    
    case 'inputdata3'
        % only read NIfTI file for mask
        [fileName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select a NIfTI file');

        if pathDir ~= 0
            % set input edit field for display
            set(h.dataIO.edit.inputData3,    'String',fullfile(pathDir,fileName));
            
            % set input directory field empty
            set(h.dataIO.edit.input,    'String',[]);

        end
        
    case 'header'
        % only read NIfTI file for mask
        [fileName,pathDir] = uigetfile({'*.mat','QSM hub header file (*.mat)'},'Select a header file');

        if pathDir ~= 0
            % set input edit field for display
            set(h.dataIO.edit.inputHeader,    'String',fullfile(pathDir,fileName));
            
            % set input directory field empty
            set(h.dataIO.edit.input,    'String',[]);

        end
end

end