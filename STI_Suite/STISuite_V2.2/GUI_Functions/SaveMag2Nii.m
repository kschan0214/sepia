function SaveMag2Nii(hObject,handles,control)

filters = {'*.nii','NII-files (*.nii)'};
switch control
    case 1
        try
            seed=get(handles.loadniiname1,'string');
        catch
            seed = 'M_.nii';
        end
        [fn,pn,filterindex] = uiputfile(filters, sprintf('Save Workspace Variables'), seed);
        STR=get(handles.LeftResolution,'string');
        eval(['resolution=' STR ';']);
        if fn~=0
            ezsave_nii(handles.Mi,strcat(pn,fn),[],resolution);
        end
    case 2
        try
            seed=get(handles.loadniiname2,'string');
        catch
            seed = 'M_.nii';
        end
        STR=get(handles.LeftResolution,'string');
        eval(['resolution=' STR ';']);
        [fn,pn,filterindex] = uiputfile(filters, sprintf('Save Workspace Variables'), seed);
        if fn~=0
            ezsave_nii(double(handles.M2),strcat(pn,fn),[],resolution);
        end
end
