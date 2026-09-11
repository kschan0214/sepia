%% relPath = get_relative_source_path(baseDir, fullPath)
%
% Input
% --------------
% baseDir       : the output directory the JSON sidecar will live in
% fullPath      : full filename of a source (input) file
%
% Output
% --------------
% relPath       : fullPath expressed relative to baseDir when possible,
%                 otherwise just the file's basename
%
% Description: Small helper to produce a tidy 'Sources' entry for a BIDS
%              JSON sidecar without depending on the source file
%              necessarily living under baseDir.
%
% Kwok-shing Chan @ DCCN
% kwokshing.chan@donders.ru.nl
% Date created: 23 August 2026 (v1.3.0)
%
function relPath = get_relative_source_path(baseDir, fullPath)

if isempty(fullPath)
    relPath = '';
    return
end

baseDir_norm  = strrep(baseDir, '\', '/');
fullPath_norm = strrep(fullPath, '\', '/');

if ~endsWith(baseDir_norm, '/')
    baseDir_norm = [baseDir_norm '/'];
end

if startsWith(fullPath_norm, baseDir_norm)
    relPath = extractAfter(fullPath_norm, baseDir_norm);
else
    [~, name, ext] = fileparts(fullPath);
    relPath = [name ext];
end

end
