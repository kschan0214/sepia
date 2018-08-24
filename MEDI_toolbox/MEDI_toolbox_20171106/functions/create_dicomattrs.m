function create_dicomattrs( )
%CREATE_DICOMATTRS Summary of this function goes here
%   Detailed explanation goes here

[STR, NAM, EXT] = fileparts(mfilename('fullpath'));
target=fullfile(STR,'dicomattrs.m');
if 2 ~= exist(target)
    warning(['Creating ' target ])
    source=fullfile(matlabroot,'toolbox','images','iptformats','dicominfo.m');
    fileID=fopen(source);T=textscan(fileID,'%s','Delimiter', '\n'); fclose(fileID);
    T1=regexprep(T{1},'^[ \t]*metadata = ', '%metadata = ');
    T2=regexprep(T1,'^[ \t]*[metadata,attrNames', ['metadata = attrs;\n%[metadata,attrNames']);
    T3=regexprep(T2,'^[ \t]*pixelDataField = ', 'pixelDataField = '''';%');
    fileID = fopen(target,'w'); fprintf(fileID,'%s\n',T3{:}); fclose(fileID);
    extrafiles={'dicom_findActualFilename.m','dicom_getFileDetails.m'};
    for j=1:length(extrafiles)
        fullf=fullfile(matlabroot,'toolbox','images','iptformats','private',extrafiles{j});
        targf=fullfile(STR,extrafiles{j});
        if 2==exist(fullf)
            warning(['Creating ' targf])
            copyfile(fullf,targf)
        end
    end
%     extraidx=zeros(size(extrafiles));
%     for i=1:length(T3)
%         for j=1:length(extrafiles)
%             if ~isempty(regexp(T3{i},extrafiles{j}))
%                 extraidx(j)=1
%             end
%         end
%     end
%     for j=find(extraidx)
%         disp(extrafiles(j))
%     end
end

target2 = fullfile(STR,['dicomparse.',mexext]);
if 3 ~= exist(target2)
    warning(['Creating ' target2 ])
    copyfile(fullfile(matlabroot,'toolbox','images','iptformats','private',['dicomparse.',mexext]),target2)
end

end

