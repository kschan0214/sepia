%% paths = generate_synthetic_phantom(outputDir)
%
% Description: deterministically builds a small multi-echo GRE-like
% synthetic phantom (magnitude + phase NIfTI, a brain mask NIfTI, and a
% sepia_header.mat) under outputDir, for use as Tier-1 regression-test
% input. Fully self-contained: no external toolbox required.
%
% Forward model: a handful of closed-form spherical susceptibility
% sources inside a synthetic "brain" ellipsoid are convolved with the
% dipole kernel (misc/qsm_algorithm/SingleStep/kspace_kernel.m) to obtain
% a local (tissue) field, a smooth polynomial background field is added,
% and multi-echo complex signal is simulated from a T2*-decay model with
% seeded complex Gaussian noise (so the resulting phase is genuinely
% wrapped and meaningfully exercises phase-unwrapping methods).
%
% Deterministic: seeded RNG (42). Safe to call repeatedly - overwrites
% any existing files of the same name under outputDir.
%
% Input
% --------------
% outputDir : directory to write the phantom files into (created if missing)
%
% Output
% --------------
% paths : struct with fields
%   .input           : 1x4 struct array matching sepiaIO's `input` argument
%                       (input(1)=phase, input(2)=magnitude, input(3)=weights
%                        [empty, unused], input(4)=sepia_header)
%   .magnitude, .phase, .mask, .header : individual file paths
%   .groundTruthChi  : .mat file with the ground-truth susceptibility map,
%                       brain mask, voxelSize and matrixSize (debug/reference
%                       use only - not consumed by any SEPIA wrapper)
%
function paths = generate_synthetic_phantom(outputDir)

if nargin < 1 || isempty(outputDir)
    error('generate_synthetic_phantom:missingOutputDir', 'An outputDir must be provided.');
end
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

rng(42, 'twister');

%% geometry / acquisition parameters
matrixSize = [32 32 24];
voxelSize  = [2 2 2];      % mm
B0         = 3;            % T
B0_dir     = [0 0 1];
gyro       = 42.57747892;  % MHz/T
nEchoes    = 5;
TE         = (3:3:3*nEchoes) * 1e-3; % s: 3,6,9,12,15 ms
delta_TE   = TE(2) - TE(1);
CF         = B0 * gyro * 1e6;

[x,y,z] = ndgrid(1:matrixSize(1), 1:matrixSize(2), 1:matrixSize(3));
cx = (matrixSize(1)+1)/2; cy = (matrixSize(2)+1)/2; cz = (matrixSize(3)+1)/2;

% synthetic "brain" ellipsoid mask
brainMask = ((x-cx)/(matrixSize(1)*0.42)).^2 + ((y-cy)/(matrixSize(2)*0.42)).^2 ...
          + ((z-cz)/(matrixSize(3)*0.42)).^2 <= 1;

%% ground-truth susceptibility sources (ppm): a few spheres of distinct sign/value
chi_true = zeros(matrixSize);
spheres = { ...
    struct('c',[cx-6, cy,   cz  ], 'r', 3.0, 'chi', +0.08), ...
    struct('c',[cx+6, cy,   cz  ], 'r', 3.0, 'chi', -0.05), ...
    struct('c',[cx,   cy-6, cz  ], 'r', 2.5, 'chi', +0.04), ...
    struct('c',[cx,   cy+6, cz+2], 'r', 2.5, 'chi', -0.03) ...
    };
for k = 1:numel(spheres)
    s = spheres{k};
    d2 = (x-s.c(1)).^2 + (y-s.c(2)).^2 + (z-s.c(3)).^2;
    chi_true(d2 <= s.r^2) = s.chi;
end
chi_true = chi_true .* brainMask;

%% forward model: dipole convolution -> local (tissue) field, in ppm
addpath(fullfile(sepiatest.sepia_home(), 'misc', 'qsm_algorithm', 'SingleStep'));
FOV = voxelSize .* matrixSize;
D   = kspace_kernel(FOV, matrixSize);           % centered (DC at index N/2+1)
fieldLocal_ppm = real(ifftn(fftn(chi_true) .* fftshift(D)));

% smooth low-order background field (dominates near/outside the mask edge,
% removable by background-field-removal methods), zero-mean inside the mask
bg = 0.30*((x-cx)/matrixSize(1)) + 0.20*((y-cy)/matrixSize(2)).^2 - 0.15*((z-cz)/matrixSize(3));
backgroundField_ppm = bg - mean(bg(brainMask(:)));

totalField_ppm = fieldLocal_ppm + backgroundField_ppm;

%% multi-echo magnitude/phase simulation with seeded complex Gaussian noise
r2 = (x-cx).^2 + (y-cy).^2 + (z-cz).^2;
S0 = 200 * exp(-r2 / (2*(matrixSize(1)*0.5)^2)) .* brainMask + 5; % small signal floor outside mask
R2s_true = 20 + 5*sin(x/3) .* brainMask;                          % 1/s, deterministic spatial variation
sigma = 3; % noise std, same units as magnitude

magnitude = zeros([matrixSize nEchoes]);
phase     = zeros([matrixSize nEchoes]);
for e = 1:nEchoes
    magn_true  = S0 .* exp(-TE(e) .* R2s_true);
    phase_true = 2*pi*gyro*1e6*B0 .* (totalField_ppm*1e-6) .* TE(e);
    noise = sigma * (randn(matrixSize) + 1i*randn(matrixSize));
    complex_signal = magn_true .* exp(1i*phase_true) + noise;
    magnitude(:,:,:,e) = abs(complex_signal);
    phase(:,:,:,e)     = angle(complex_signal); % naturally wrapped to (-pi,pi]
end

%% write NIfTI outputs
% NOTE: save_nii_quick.m (used elsewhere in SEPIA) requires an "untouch"
% style template obtained from load_untouch_nii(existingFile) - there is
% no pre-existing file to load here, so this from-scratch generator uses
% the NIfTI toolbox's own make_nii/save_nii pair instead (bundled at
% utils/nifti/NIfTI_20140122/, already on path via sepia_addpath).
paths = struct();
paths.magnitude      = fullfile(outputDir, 'phantom_mag.nii.gz');
paths.phase          = fullfile(outputDir, 'phantom_phase.nii.gz');
paths.mask           = fullfile(outputDir, 'phantom_mask.nii.gz');
paths.header         = fullfile(outputDir, 'phantom_header.mat');
paths.groundTruthChi = fullfile(outputDir, 'phantom_groundtruth_chi.mat');

save_nii(make_nii(single(magnitude), voxelSize, [], 16), paths.magnitude);
save_nii(make_nii(single(phase),     voxelSize, [], 16), paths.phase);
save_nii(make_nii(int16(brainMask),  voxelSize, [], 4),  paths.mask); % datatype 4 = int16, binary mask

%% write sepia_header.mat directly (same 7 loose variables save_sepia_header.m
% itself writes - see utils/save_sepia_header.m's final save() call - since
% that function otherwise expects to read an existing NIfTI/DICOM/JSON
% dataset on disk rather than pure programmatic values)
save(paths.header, 'voxelSize', 'matrixSize', 'CF', 'delta_TE', 'TE', 'B0_dir', 'B0');

save(paths.groundTruthChi, 'chi_true', 'brainMask', 'voxelSize', 'matrixSize');

%% assemble sepiaIO-style input struct array
input(1).name = paths.phase;
input(2).name = paths.magnitude;
input(3).name = '';
input(4).name = paths.header;
paths.input = input;

end
