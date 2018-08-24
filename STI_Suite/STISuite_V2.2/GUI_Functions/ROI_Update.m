function ROI_Update(hObject,handles)


TempROIvalues=ones(size(handles.M2));
myvarname=get(handles.RightImageSave,'string');
myvarname=[myvarname get(handles.CoilNum,'string')];
MyROIAnalysis=evalin('base', myvarname);
MyROIs=MyROIAnalysis.ROI;

for ii=1:9
    TempROIvalues(MyROIs{ii})=ii*0.01+handles.baseAlpha;
end            
handles.M2=TempROIvalues;

guidata(gcbo,handles);
