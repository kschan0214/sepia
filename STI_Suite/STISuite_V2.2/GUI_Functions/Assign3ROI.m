function Assign3ROI(hObject,handles,ctrl)
switch ctrl 
    case 1
            TempROIvalues=ones(size(handles.M2));
            myvarname=get(handles.edit17,'string');
            MyROIAnalysis=evalin('base', myvarname);
            for ii=1:9
                MyROIs=MyROIAnalysis==ii;
                TempROIvalues(MyROIs)=ii*0.01+handles.baseAlpha;
            end            
            handles.M2=TempROIvalues;
            guidata(gcbo,handles);
    
end