%% save_json_sidecar(niiFilename, info)
%
% Input
% --------------
% niiFilename   : full filename of the NIfTI output the sidecar describes
%                 ('*.nii' or '*.nii.gz'). The sidecar is written next to
%                 it with the same basename and a '.json' extension
% info          : (optional) structure with the extra fields to write
%                 into the sidecar, e.g.
%                   info.Description = 'Local field map after V-SHARP';
%                   info.Units       = 'Hz';
%                   info.Sources     = {'sub-01_part-phase_desc-unwrapped.nii.gz'};
%                   info.Parameters  = algorParam.bfr;
%                 'Units' should always be set for quantitative maps
%                 (e.g. 'Hz', 'rad', 'ppm', '1/s', 's', 'arbitrary') and
%                 omitted for masks/binary images.
%
% Description: Writes a BIDS-Derivatives-style JSON sidecar next to a
%              SEPIA NIfTI output for provenance/BIDS compatibility.
%              Common fields (SoftwareName, SoftwareVersion,
%              GeneratedDate) are always added automatically.
%
% Kwok-shing Chan @ DCCN
% kwokshing.chan@donders.ru.nl
% Date created: 23 August 2026 (v1.3.0)
%
function save_json_sidecar(niiFilename, info)

if nargin < 2 || isempty(info)
    info = struct();
end

sepia_universal_variables;

% common provenance fields go first
jsonStruct = struct();
jsonStruct.SoftwareName    = 'SEPIA';
jsonStruct.SoftwareVersion = SEPIA_version;
jsonStruct.GeneratedDate   = datestr(datetime('now'),'yyyy-mm-ddTHH:MM:SS');

% append/overwrite with caller-supplied fields (e.g. Description, Units,
% Sources, Parameters)
fields = fieldnames(info);
for k = 1:numel(fields)
    jsonStruct.(fields{k}) = info.(fields{k});
end

jsonFilename = get_json_filename_from_nifti(niiFilename);

try
    txt = jsonencode(jsonStruct, 'PrettyPrint', true);
catch
    % 'PrettyPrint' option requires MATLAB R2021a+; fall back gracefully
    txt = jsonencode(jsonStruct);
end

fid = fopen(jsonFilename,'w');
if fid == -1
    warning('save_json_sidecar:cannotWrite','Cannot open ''%s'' for writing. JSON sidecar was not saved.', jsonFilename);
    return
end
fwrite(fid, txt);
fclose(fid);

end

%% derive '<basename>.json' from a '.nii'/'.nii.gz' filename
function jsonFilename = get_json_filename_from_nifti(niiFilename)

if endsWith(niiFilename, '.nii.gz')
    jsonFilename = [niiFilename(1:end-numel('.nii.gz')) '.json'];
elseif endsWith(niiFilename, '.nii')
    jsonFilename = [niiFilename(1:end-numel('.nii')) '.json'];
else
    [pathstr, name] = fileparts(niiFilename);
    jsonFilename = fullfile(pathstr, [name '.json']);
end

end
