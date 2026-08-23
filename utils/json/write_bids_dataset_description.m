%% write_bids_dataset_description(outputDir)
%
% Input
% --------------
% outputDir     : output directory of a SEPIA processing pipeline
%
% Description: Writes a minimal BIDS-Derivatives 'dataset_description.json'
%              at the root of a SEPIA output directory, if one does not
%              already exist. Existing files are left untouched so this
%              is safe to call every time an output pipeline runs.
%
% Kwok-shing Chan @ DCCN
% kwokshing.chan@donders.ru.nl
% Date created: 23 August 2026 (v1.3.0)
%
function write_bids_dataset_description(outputDir)

descriptionFilename = fullfile(outputDir, 'dataset_description.json');

% do not overwrite a description that already exists (e.g. user-edited,
% or from a previous run)
if exist(descriptionFilename, 'file') == 2
    return
end

sepia_universal_variables;

description                 = struct();
description.Name            = 'SEPIA outputs';
description.BIDSVersion     = '1.9.0';
description.DatasetType     = 'derivative';
description.GeneratedBy     = {struct('Name','SEPIA','Version',SEPIA_version,'CodeURL','https://github.com/kschan0214/sepia')};

try
    txt = jsonencode(description, 'PrettyPrint', true);
catch
    txt = jsonencode(description);
end

fid = fopen(descriptionFilename,'w');
if fid == -1
    warning('write_bids_dataset_description:cannotWrite','Cannot open ''%s'' for writing.', descriptionFilename);
    return
end
fwrite(fid, txt);
fclose(fid);

end
