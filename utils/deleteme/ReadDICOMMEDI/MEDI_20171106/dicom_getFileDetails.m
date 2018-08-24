function details = dicom_getFileDetails(filename)
%dicom_getFileDetails   Get file's full pathname, size, etc.

%   Copyright 2016 The MathWorks, Inc.

% Get the fully qualified path to the file.  
fid = fopen(filename);
if (fid < 0)
    augmentedFilename = dicom_findActualFilename(filename);
    
    if (~isempty(augmentedFilename))
        fid = fopen(augmentedFilename);
    else
        error(message('images:dicomread:fileNotFound', filename))
    end
end

fullPathFilename = fopen(fid);
details.name = fullPathFilename;

details.bytes = getFileLength(fid);

fclose(fid);

d = dir(fullPathFilename);
details.date = d.date;

end


function fileLength = getFileLength(fid)

fseek(fid, 0, 'eof');
fileLength = ftell(fid);

end