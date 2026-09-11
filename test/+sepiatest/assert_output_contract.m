%% sepiatest.assert_output_contract(testCase, outputPrefix, suffix, expectations, maskFile)
%
% Description: asserts a set of non-numeric regression properties about
% SEPIA's output for one pipeline run:
%   1. each expected output NIfTI file exists
%   2. its BIDS-Derivatives-style JSON sidecar is present/absent exactly
%      as expected (see save_json_sidecar.m call sites)
%   3. no NaN/Inf voxel inside maskFile (if provided)
%   4. matrix size matches the mask's matrix size (if provided)
%
% Input
% --------------
% testCase     : a matlab.unittest.qualifications.Qualifiable
% outputPrefix : the same '<dir>/<prefix>' string passed to sepiaIO as `output`
%                (filenames are '<dir>/<prefix>_<key><suffix>', matching
%                 wrapper/private/construct_output_filename.m's convention)
% suffix       : nifti extension, e.g. '.nii.gz' (see get_nifti_extension_from_input.m)
% expectations : cell array of {key, expectJson} rows, e.g.
%                {{'Chimap', true}, {'localfield', true}, {'fieldmap', true}}
% maskFile     : (optional) mask NIfTI filename to check matrix size / NaN-Inf against;
%                pass '' to skip that check
%
function assert_output_contract(testCase, outputPrefix, suffix, expectations, maskFile)

if nargin < 5
    maskFile = '';
end

[outDir, prefixName] = fileparts(outputPrefix);
% fileparts on '/dir/sepia' with no extension gives prefixName='sepia';
% construct_output_filename appends '_' + key + suffix to '<prefixName>_'
baseName = fullfile(outDir, [prefixName '_']);

maskImg = [];
if ~isempty(maskFile)
    testCase.assertTrue(isfile(maskFile), sprintf('mask file missing: %s', maskFile));
    maskImg = load_nii_img_only(maskFile) > 0;
end

for k = 1:numel(expectations)
    key = expectations{k}{1};
    expectJson = expectations{k}{2};

    niiFile = [baseName key suffix];
    testCase.verifyTrue(isfile(niiFile), sprintf('Expected output file missing: %s', niiFile));
    if ~isfile(niiFile)
        continue % nothing further to check for this key
    end

    jsonFile = json_sidecar_name(niiFile);
    if expectJson
        testCase.verifyTrue(isfile(jsonFile), sprintf('Expected JSON sidecar missing: %s', jsonFile));
    else
        testCase.verifyFalse(isfile(jsonFile), sprintf('Unexpected JSON sidecar present (not documented): %s', jsonFile));
    end

    img = load_nii_img_only(niiFile);
    if ~isempty(maskImg)
        if isequal(size(img,1:3), size(maskImg,1:3))
            v = img(repmat(maskImg,[1 1 1 size(img,4)]));
            testCase.verifyFalse(any(isnan(v)), sprintf('%s: NaN voxel(s) found inside mask', key));
            testCase.verifyFalse(any(isinf(v)), sprintf('%s: Inf voxel(s) found inside mask', key));
        else
            testCase.verifyEqual(size(img,1:3), size(maskImg,1:3), ...
                sprintf('%s: output matrix size does not match input mask matrix size', key));
        end
    end
end

end

%% mirror get_json_filename_from_nifti's logic in save_json_sidecar.m
function jsonFilename = json_sidecar_name(niiFilename)

if endsWith(niiFilename, '.nii.gz')
    jsonFilename = [niiFilename(1:end-numel('.nii.gz')) '.json'];
elseif endsWith(niiFilename, '.nii')
    jsonFilename = [niiFilename(1:end-numel('.nii')) '.json'];
else
    [pathstr, name] = fileparts(niiFilename);
    jsonFilename = fullfile(pathstr, [name '.json']);
end

end
