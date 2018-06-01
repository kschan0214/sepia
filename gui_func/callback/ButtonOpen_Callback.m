function ButtonOpen_Callback(source,eventdata,h,field)
% get directory and display it on GUI

% global h

% output base name
prefix = 'squirrel';

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
        end
        
    case 'output'
        
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.dataIO.edit.output,'String',[pathDir filesep prefix]);
        end
end

end