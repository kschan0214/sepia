function metadata = dicominfo(filename, varargin)
%DICOMINFO  Read metadata from DICOM message.
%   INFO = DICOMINFO(FILENAME) reads the metadata from the compliant
%   DICOM file specified in the string FILENAME.
%
%   INFO = DICOMINFO(FILENAME, 'dictionary', D) uses the data dictionary
%   file given in the string D to read the DICOM message.  The file in D
%   must be on the MATLAB search path.  The default value is dicom-dict.mat.
%
%   INFO = DICOMINFO(..., 'UseVRHeuristic', TF) instructs the parser to
%   use a heuristic to help read certain noncompliant files which switch
%   value representation (VR) modes incorrectly. A warning will be
%   displayed if the heuristic is employed. When TF is true (the default),
%   a small number of compliant files will not be read correctly. Set TF to
%   false to read these compliant files.
%
%   INFO = DICOMINFO(..., 'UseDictionaryVR', TF) specifies whether the
%   datatypes in INFO should conform to the data dictionary, regardless of
%   what information is present in the file. The default value of TF is
%   false, which uses the file's VR codes even if they differ from the data
%   dictionary. Most of the time it is unnecessary to set the field, since
%   file contents and the data dictionary almost always agree. When TF is
%   false (the default) a warning is issued when they do not agree. Set TF
%   to true when the warning is issued and providing INFO to DICOMWRITE
%   causes errors.
%
%   Example:
%
%     info = dicominfo('CT-MONO2-16-ankle.dcm');
%
%   See also DICOMREAD, DICOMWRITE, DICOMDISP, DICOMDICT, DICOMUID.

%   Copyright 1993-2016 The MathWorks, Inc.


% Parse input arguments.
if (nargin < 1)

error(message('images:dicominfo:tooFewInputs'))

end

% Set the dictionary.
args = parseInputs(varargin{:});
dicomdict('set_current', args.Dictionary)
dictionary = dicomdict('get_current');

% Get the metadata.
try

fileDetails = dicom_getFileDetails(filename);

% Ensure the file is actually DICOM.
if (~isdicom(fileDetails.name))
error(message('images:dicominfo:notDICOM'))
end

% Parse the DICOM file.
attrs = dicomparse(fileDetails.name, ...
fileDetails.bytes, ...
getMachineEndian, ...
false, ...
dictionary, ...
args.UseVRHeuristic); 

% Process the raw attributes.
metadata = attrs;
%[metadata,attrNames] = processMetadata(attrs, true, dictionary, args.UseDictionaryVR);
%metadata = dicom_set_imfinfo_values(metadata);
%metadata = setMoreImfinfoValues(metadata, fileDetails);
%metadata = processOverlays(metadata, dictionary);
%metadata = processCurveData(metadata, attrs, attrNames);

pixelDataField = '';%dicomlookup_actions('7fe0', '0010', dictionary);
if (isfield(metadata, pixelDataField))
%metadata = rmfield(metadata, pixelDataField);
end

catch ME

dicomdict('reset_current');
rethrow(ME)

end

% Reset the dictionary.
dicomdict('reset_current');



function [metadata,attrNames] = processMetadata(attrs, isTopLevel, dictionary, useDictionaryVR)

if (isempty(attrs))
%metadata = [];
return
end

% Create a structure for the output and get the names of attributes.
[metadata, attrNames] = createMetadataStruct(attrs, isTopLevel, dictionary);

% Fill the metadata structure, converting data along the way.
for currentAttr = 1:numel(attrNames)

this = attrs(currentAttr);
metadata.(attrNames{currentAttr}) = convertRawAttr(this, dictionary, useDictionaryVR, metadata);

end



function processedAttr = convertRawAttr(rawAttr, dictionary, useDictionaryVR, metadata)

% Information about whether to swap is contained in the attribute.
swap = needToSwap(rawAttr);

[dictionaryVR, name] = findVRFromTag(rawAttr.Group, rawAttr.Element, dictionary);

% Determine the correct output encoding.
if (isempty(rawAttr.VR))

% Look up VR for implicit VR files.  Use 'UN' for unknown
% tags.  (See PS 3.5 Sec. 6.2.2.)
if (~isempty(dictionaryVR))

% Some attributes have a conditional VR.  Pick the first.
rawAttr.VR = dictionaryVR;
if (numel(rawAttr.VR) > 2)
rawAttr.VR = rawAttr.VR(1:2);
end

else
rawAttr.VR = 'UN';
end

else

if (~isempty(dictionaryVR) && ...
~isequal(dictionaryVR, 'UN') && ...
isempty(strfind(dictionaryVR, rawAttr.VR)) && ...
(rem(rawAttr.Group, 2) ~= 1))

if (~useDictionaryVR)
attrTag = sprintf('(%04X,%04X)', rawAttr.Group, rawAttr.Element);
warning(message('images:dicominfo:fileVRDoesNotMatchDictionary', ...
name, attrTag, dictionaryVR, rawAttr.VR))
else
rawAttr.VR = dictionaryVR;
end

end

end

% Convert raw data.  (See PS 3.5 Sec. 6.2 for full VR details.)
switch (rawAttr.VR)
case  {'AE', 'AS', 'CS', 'DA', 'DT', 'TM', 'UI', 'UR'}

processedAttr = deblankAndStripNulls(char(rawAttr.Data));

case {'AT'}

% For historical reasons don't transpose AT.
processedAttr = dicom_typecast(rawAttr.Data, 'uint16', swap);

case {'DS', 'IS'}

processedAttr = sscanf(char(rawAttr.Data), '%f\\');

case {'FL', 'OF'}

processedAttr = dicom_typecast(rawAttr.Data, 'single', swap)';

case {'FD', 'OD'}

processedAttr = dicom_typecast(rawAttr.Data, 'double', swap)';

case {'LO', 'LT', 'SH', 'ST', 'UC', 'UT'}

specificCharacterSet = dicom_get_SpecificCharacterSet(metadata, dictionary);
processedAttr = getUnicodeStringFromBytes(rawAttr.Data, specificCharacterSet{end});
processedAttr = deblankAndStripNulls(processedAttr);

case 'OB'

processedAttr = rawAttr.Data';

case {'OW', 'US'}

processedAttr = dicom_typecast(rawAttr.Data, 'uint16', swap)';

case 'PN'

processedAttr = parsePerson(rawAttr.Data, metadata, dictionary);

case 'SL'

processedAttr = dicom_typecast(rawAttr.Data, 'int32', swap)';

case 'SQ'

processedAttr = parseSequence(rawAttr.Data, dictionary, useDictionaryVR);

case 'SS'

processedAttr = dicom_typecast(rawAttr.Data, 'int16', swap)';

case 'UL'

processedAttr = dicom_typecast(rawAttr.Data, 'uint32', swap)';

case 'UN'

% It's possible that the attribute contains a private sequence
% with implicit VR; in which case the Data field contains the
% parsed sequence.
if (isstruct(rawAttr.Data))
processedAttr = parseSequence(rawAttr.Data, dictionary, useDictionaryVR);
else
processedAttr = rawAttr.Data';
end

otherwise

% PS 3.5-1999 Sec. 6.2 indicates that all unknown VRs can be
% interpretted as UN.  
processedAttr = rawAttr.Data';

end

% Change empty arrays to 0-by-0.
if isempty(processedAttr)
processedAttr = reshape(processedAttr, [0 0]);
end



function byteOrder = getMachineEndian

persistent endian

if (~isempty(endian))
byteOrder = endian;
return
end

[~, ~, endian] = computer;
byteOrder = endian;



function args = parseInputs(varargin)

% Set default values
args.Dictionary = dicomdict('get');
args.UseVRHeuristic = true;
args.UseDictionaryVR = false;

% Parse arguments based on their number.
if (nargin > 0)

paramStrings = {'dictionary', 'usevrheuristic', 'usedictionaryvr'};

% For each pair
for k = 1:2:length(varargin)
param = lower(varargin{k});

if (~ischar(param))
error(message('images:dicominfo:parameterNameNotString'));
end

idx = dicom_strmatch(param, paramStrings);

if (isempty(idx))
error(message('images:dicominfo:unrecognizedParameterName', param));
elseif (length(idx) > 1)
error(message('images:dicominfo:ambiguousParameterName', param));
end

switch (paramStrings{idx})
case 'dictionary'

if (k == length(varargin))
error(message('images:dicominfo:missingDictionary'));
else
args.Dictionary = varargin{k + 1};
end

case 'usevrheuristic'

if (k == length(varargin))
error(message('images:dicominfo:missingVRHeuristic'));
else
args.UseVRHeuristic = varargin{k+1};
validateattributes(args.UseVRHeuristic, {'logical'}, ...
{'scalar', 'nonempty'}, mfilename, 'UseVRHeuristic', k)
end

case 'usedictionaryvr'

if (k == length(varargin))
error(message('images:dicominfo:missingUseDictionaryVR'));
else
args.UseDictionaryVR = varargin{k+1};
validateattributes(args.UseDictionaryVR, {'logical'}, ...
{'scalar', 'nonempty'}, mfilename, 'UseDictionaryVR', k)
end

end

end

end



function personName = parsePerson(rawData, metadata, dictionary)
%PARSEPERSON  Get the various parts of a person name

% A description and examples of PN values is in PS 3.5-2000 Table 6.2-1.

pnParts = {'FamilyName'
'GivenName'
'MiddleName'
'NamePrefix'
'NameSuffix'};

if (isempty(rawData))
personName = makePerson(pnParts);
return
end

splitRawData = tokenizeRawData(rawData, '\');

specificCharacterSet = dicom_get_SpecificCharacterSet(metadata, dictionary);

personName = struct([]);

for p = 1:numel(splitRawData)
% Extend the structure for this part of the PN attribute.
personName(p).(pnParts{1}) = '';

if (hasMultiplePNComponents(splitRawData{p}))
unicodeString = convertMultipleComponentsToUnicode(splitRawData{p}, ...
specificCharacterSet);
else
unicodeString = getUnicodeStringFromBytes(splitRawData{p}, ...
specificCharacterSet{1});
end
unicodeString = deblankAndStripNulls(unicodeString);

% Alphabetic, ideographic, and phonetic components are separated by '='.
components = tokenize(unicodeString, '=');

if (isempty(components))
continue
end

for componentIdx = 1:numel(components)
% Get the separate parts of the person's name from the component.
componentParts = tokenize(components{componentIdx}, '^');

% The DICOM standard requires that PN values have five or fewer
% values separated by "^".  Some vendors produce files with more
% than these person name parts.
if (numel(componentParts) <= 5)
% If there are the correct numbers, put them in separate fields.
for q = 1:length(componentParts)

if (isfield(personName, pnParts{q}) && ~isempty(personName(p).(pnParts{q})))
personName(p).(pnParts{q}) = [personName(p).(pnParts{q}) '=' componentParts{q}];
else
personName(p).(pnParts{q}) = componentParts{q};
end
end
else
% If there are more, just return the whole string.
personName(p).FamilyName = unicodeString;
break
end
end
end

for i = 1:numel(personName)
for j = 1:numel(fieldnames(personName))
if (isempty(personName(i).(pnParts{j})))
personName(i).(pnParts{j}) = '';
end
end
end


function splitRawData = tokenizeRawData(rawData, delimiter)

assert(~isempty(rawData))

delimiterLocations = find(rawData == uint8(delimiter));

if isempty(delimiterLocations)
splitRawData = {rawData};
return
end

splitRawData = cell(1, numel(delimiterLocations) + 1);

start = 1;
for p = 1:numel(delimiterLocations)
stop = delimiterLocations(p) - 1;
if (stop >= start)
splitRawData{p} = rawData(start:stop);
end
start = stop + 2;
end

splitRawData{end} = rawData(start:end);


function TF = hasMultiplePNComponents(rawData)

if (~isempty(rawData))
TF = ~isempty(find(rawData == '=', 1));
else
TF = false;
end


function personStruct = makePerson(pnParts)
%MAKEPERSON  Make an empty struct containing the PN fields.

for p = 1:numel(pnParts)
personStruct.(pnParts{p}) = '';
end


function unicodeString = convertMultipleComponentsToUnicode(rawData, specificCharacterSet)

% Components are divided by equal signs ('=')
separatorIndices = find(rawData == uint8('='));

unicodeString = '';

start = 1;
for p = 1:numel(separatorIndices)

if (p == 1)
targetCharSet = specificCharacterSet{1};
else
targetCharSet = specificCharacterSet{end};
end

stop = separatorIndices(p) - 1;
substring = rawData(start:stop);
substring = getUnicodeStringFromBytes(substring, targetCharSet);
unicodeString = [unicodeString '=' substring]; %#ok<AGROW>

start = separatorIndices(p) + 1;
end

substring = rawData(start:end);
substring = getUnicodeStringFromBytes(substring, specificCharacterSet{end});
unicodeString = [unicodeString '=' substring];
unicodeString(1) = '';  % Remove leading '='


function processedStruct = parseSequence(attrs, dictionary, useDictionaryVR)

numItems = countItems(attrs);
itemNames = getItemNames(numItems);

% Initialize the structure to contain this structure.
structInitializer = cat(1, itemNames, cell(1, numItems));
processedStruct = struct(structInitializer{:});

% Process each item (but not delimiters).
item = 0;
for idx = 1:numel(attrs)

this = attrs(idx);
if (~isDelimiter(this))
item = item + 1;
processedStruct.(itemNames{item}) = processMetadata(this.Data, false, dictionary, useDictionaryVR);
end

end



function header = getImfinfoFields

header = {'Filename',      ''
'FileModDate',   ''
'FileSize',      []
'Format',        'DICOM'
'FormatVersion', 3.0
'Width',         []
'Height',        []
'BitDepth',      []
'ColorType',     ''}';



function metadata = setMoreImfinfoValues(metadata, d)

metadata.Filename    = d.name;
metadata.FileModDate = d.date;
metadata.FileSize    = d.bytes;



function [metadata, attrNames] = createMetadataStruct(attrs, isTopLevel, dictionary)

% Get the attribute names.
totalAttrs = numel(attrs);
attrNames = cell(1, totalAttrs);

for currentAttr = 1:totalAttrs
attrNames{currentAttr} = ...
dicomlookup_actions(attrs(currentAttr).Group, ...
attrs(currentAttr).Element, ...
dictionary);

% Empty attributes indicate that a public/retired attribute was
% not found in the data dictionary.  This used to be an error
% condition, but is easily resolved by providing a special
% attribute name.
if (isempty(attrNames{currentAttr}))
attrNames{currentAttr} = sprintf('Unknown_%04X_%04X', ...
attrs(currentAttr).Group, ...
attrs(currentAttr).Element);
end
end

% Remove duplicate attribute names.  Keep the last appearance of the attribute.
[tmp, reorderIdx] = unique(attrNames);
if (numel(tmp) ~= totalAttrs)
warning(message('images:dicominfo:attrWithSameName'))
end

uniqueAttrNames = attrNames(sort(reorderIdx));
uniqueTotalAttrs = numel(uniqueAttrNames);

% Create a metadata structure to hold the parsed attributes.  Use a
% cell array initializer, which has a populated section for IMFINFO
% data and an unitialized section for the attributes from the DICOM
% file.
if (isTopLevel)
structInitializer = cat(2, getImfinfoFields(), ...
cat(1, uniqueAttrNames, cell(1, uniqueTotalAttrs)));
else
structInitializer = cat(1, uniqueAttrNames, cell(1, uniqueTotalAttrs));
end

%metadata = struct(structInitializer{:});



function str = deblankAndStripNulls(str)
%DEBLANKANDDENULL  Deblank a string, treating char(0) as a blank.

if (isempty(str))
return
end

while (~isempty(str) && (str(end) == 0))
str(end) = '';
end

str = deblank(str);



function [vr, name] = findVRFromTag(group, element, dictionary)

% Look up the attribute.
attr = dicomlookup_helper(group, element, dictionary);

% Get the vr.
if (~isempty(attr))

vr = attr.VR;
name = attr.Name;

else

% Private creator attributes should be treated as CS.
if ((rem(group, 2) == 1) && (element == 0))
vr = 'UL';
elseif ((rem(group, 2) == 1) && (element < 256))
vr = 'CS';
else
vr = 'UN';
end

name = '';

end



function out = processOverlays(in, dictionary)

out = in;

% Look for overlays.
allFields = fieldnames(in);
idx = strmatch('OverlayData', allFields);

if (isempty(idx))
return
end

% Convert each overlay data attribute.
for p = 1:numel(idx)

olName = allFields{idx(p)};

% The overlay fields can be present but empty.
if (isempty(in.(olName)))
continue;
end

% Which repeating group is this?
[group, element] = dicomlookup_actions(olName, dictionary);

% Get relevant details.  All overlays are in groups 6000 - 60FE.
overlay.Rows    = double(in.(dicomlookup_actions(group, '0010', dictionary)));
overlay.Columns = double(in.(dicomlookup_actions(group, '0011', dictionary)));

sppTag = dicomlookup_actions(group, '0012', dictionary);
if (isfield(in, sppTag))
overlay.SamplesPerPixel = double(in.(sppTag));
else
overlay.SamplesPerPixel = 1;
end

bitsTag = dicomlookup_actions(group, '0100', dictionary);
if (isfield(in, bitsTag))
overlay.BitsAllocated = double(in.(bitsTag));
else
overlay.BitsAllocated = 1;
end

numTag = dicomlookup_actions(group, '0015', dictionary);
if (isfield(in, numTag))
overlay.NumberOfFrames = double(in.(numTag));
else
overlay.NumberOfFrames = 1;
end

% We could potential support more overlays later.
if ((overlay.BitsAllocated > 1) || (overlay.SamplesPerPixel > 1))

warning(message('images:dicominfo:unsupportedOverlay', sprintf( '(%04X,%04X)', group, element )));
continue;

end

% Process the overlay.
for frame = 1:(overlay.NumberOfFrames)

overlayData = tobits(in.(olName));
numSamples = overlay.Columns * overlay.Rows * overlay.NumberOfFrames;
out.(olName) = permute(reshape(overlayData(1:numSamples), ...
overlay.Columns, ...
overlay.Rows, ...
overlay.NumberOfFrames), ...
[2 1 3]);

end

end



function out = processCurveData(in, attrs, attrNames)
% Reference - PS 3.3 - 2003 C 10.2.  The Curve Data's final data type
% depends on DataValueRepresentation. Process those attributes here.

% Passing in attrs because we may need to swap the data depending on the
% endianness of the machine and the attribute.

% default
out = in;

% Look for Curve Data.
allFields = fieldnames(in);
idx = strmatch('CurveData', allFields);

if (isempty(idx))
return
end

% All the Curve Data attributes will have the same endianness, so we just need
% to check one attribute and apply that setting to the rest.

curveDataName = allFields{idx(1)};
curveDataLoc = strncmp(curveDataName,attrNames, length(curveDataName));
swap = needToSwap(attrs(curveDataLoc));

for p = 1 : numel(idx)

curveDataName = allFields{idx(p)};

underscore = strfind(curveDataName,'_');
% The data type of the Curve Data comes from the
% DataValueRepresentation attribute.
numOfRepeatedAttr = curveDataName(underscore+1:end);
dvrName = strcat('DataValueRepresentation','_',numOfRepeatedAttr);

if ~isfield(in,dvrName)
% do nothing
continue;
else

dataType = in.(dvrName);
% See PS 3.3-2003, C 10.2.1.2
switch dataType
case 0
expDataType = 'uint16';
case 1
expDataType = 'int16';
case 2
expDataType = 'single';
case 3
expDataType = 'double';
case 4
expDataType = 'int32';
otherwise
warning(message('images:dicominfo:unknownDataType', dataType));
end

% We need to undo any previous swapping before calling typecast, and do the
% final swapping afterwards.

if swap
out.(curveDataName) = ...
swapbytes(typecast(swapbytes(in.(curveDataName)), ...
expDataType));
else
out.(curveDataName) = typecast(in.(curveDataName), expDataType);
end
end
end



function itemNames = getItemNames(numberOfItems)

% Create a cell array of item names, which can be quickly used.
persistent namesCell
if (isempty(namesCell))
namesCell = generateItemNames(50);
end

% If the number of cached names is too small, expand it and recache.
if (numberOfItems > numel(namesCell))
namesCell = generateItemNames(numberOfItems);
end

% Return the first n item names.
itemNames = namesCell(1:numberOfItems);



function namesCell = generateItemNames(numberOfItems)

namesCell = cell(1, numberOfItems);
for idx = 1:numberOfItems
namesCell{idx} = sprintf('Item_%d', idx);
end



function tf = needToSwap(currentAttr)

switch (getMachineEndian)
case 'L'
if (currentAttr.IsLittleEndian)
tf = false;
else
tf = true;
end

case 'B'
if (currentAttr.IsLittleEndian)
tf = true;
else
tf = false;
end

otherwise
error(message('images:dicominfo:unknownEndian', getMachineEndian))

end



function tf = isDelimiter(attr)

% True if (FFFE,E00D) or (FFFE,E0DD).
tf = (attr.Group == 65534) && ...
((attr.Element == 57357) || (attr.Element == 57565));



function count = countItems(attrs)

if (isempty(attrs))
count = 0;
else
% Find the items (FFFE,E000) in the array of attributes (all of
% which are item tags or delimiters; no normal attributes
% appear in attrs here). 
idx = find(([attrs(:).Group] == 65534) & ...
([attrs(:).Element] == 57344));
count = numel(idx);
end


function unicodeString = getUnicodeStringFromBytes(rawData, originalCharacterSet)

icuCharacterSet = dicom_getConverterString(originalCharacterSet);

% Remove a particular control code sequence ICU doesn't handle well.
rawData = uint8(strrep(char(rawData), char([27 36 41 67]), ''));
unicodeString = native2unicode(rawData, icuCharacterSet);
