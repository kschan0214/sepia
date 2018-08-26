function Assign2BaseROI(hObject,handles,ctrl)

switch ctrl 
    case 1
            TempROIvalues=ones(size(handles.M2));
            myvarname=get(handles.RightImageLoad,'string');
            MyROIAnalysis=evalin('base', myvarname);
            MyROIs=MyROIAnalysis.ROI;
            
            for ii=1:9
                TempROIvalues(MyROIs{ii})=ii*0.01+handles.baseAlpha;
            end            
            handles.M2=TempROIvalues;
            guidata(gcbo,handles);
    case 2
        
        DataSetNumber=round(str2double(get(handles.CoilNum,'String')));
        Imag=handles.Mi(:,:,:,DataSetNumber);
        MyROIs=handles.M2;
        for ii=1:9
            ROI_ii=find(MyROIs==ii*0.01+handles.baseAlpha);
            ROI_Value_ii=Imag(ROI_ii);
            MyROIAnalysis.ROI{ii}=ROI_ii;
            MyROIAnalysis.data{ii}=ROI_Value_ii;
        end
         myvarname=get(handles.RightImageSave,'string');
        assignin('base',myvarname,MyROIAnalysis);
        
        set(handles.status,'string',[myvarname ' saved.'])
        
    case 3
        [file,file_path] = uigetfile('*.mat','MultiSelect','on');
        my_var=load([file_path file]);
        my_name=fieldnames(my_var);

        
            TempROIvalues=ones(size(handles.M2));

            MyROIAnalysis=getfield(my_var,my_name{1});
            MyROIs=MyROIAnalysis.ROI;
            
            for ii=1:9
                TempROIvalues(MyROIs{ii})=ii*0.01+handles.baseAlpha;
            end            
            handles.M2=TempROIvalues;
            guidata(gcbo,handles);

end