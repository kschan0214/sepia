function fullFilename = dicom_findActualFilename(originalFilename)
%dicom_findActualFilename  Find DICOM file's name with common extensions.

%   Copyright 2016 The MathWorks, Inc.

if (exist(originalFilename, 'file'))
    
    fullFilename = originalFilename;
    
elseif (exist([originalFilename '.dcm'], 'file'))
    
    fullFilename = [originalFilename '.dcm'];
    
elseif (exist([originalFilename '.dic'], 'file'))
    
    fullFilename = [originalFilename '.dic'];
    
elseif (exist([originalFilename '.dicom'], 'file'))
    
    fullFilename = [originalFilename '.dicom'];
    
elseif (exist([originalFilename '.img'], 'file'))
    
    fullFilename = [originalFilename '.img'];
    
else
    
    fullFilename = '';
    return
    
end
    
end
