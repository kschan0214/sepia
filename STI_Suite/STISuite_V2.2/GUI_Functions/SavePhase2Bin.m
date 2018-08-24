function SavePhase2Bin(hObject,handles)

filters = {'*.bin','BIN-files (*.bin)'};
seed = 'InputPhaseImage.bin';
[fn,pn,filterindex] = uiputfile(filters, sprintf('Save Workspace Variables'), seed);
datatype='float';
fid = fopen(strcat(pn,fn),'w','l');
fwrite(fid, handles.Mi,datatype); 
fclose(fid);
disp(size(handles.Mi))
disp('done')