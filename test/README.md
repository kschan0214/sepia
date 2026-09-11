# SEPIA regression test suite

Three tiers:

- **Tier 1 — smoke test** (`test/tier1_smoke/`): runs the full one-stop pipeline
  (`sepiaIO` → `SepiaIOWrapper`) on a small, deterministically-generated synthetic
  phantom, using only methods that need no external toolbox (phase unwrap `'None'`,
  background field removal `'VSHARP'` (the built-in variant), and QSM `'TKD'` /
  `'Closed-form solution'` / `'iLSQR'`). Runs in well under a minute. This is the
  tier that runs in CI (`.github/workflows/tier1-smoke.yml`) on every push/PR.
  Also includes `TestOddMatrixSize.m`, which runs the same toolbox-free pipeline at
  several odd/even matrix sizes to check SEPIA's zero-pad/crop handling
  (`utils/zeropad_odd_dimension.m`).
- **Tier 2 — synthetic-phantom regression matrix, all methods** (`test/tier2_matrix/`):
  one representative pipeline run per supported method **including toolbox-dependent
  ones**, run against the same small synthetic phantom Tier 1 uses (not a real
  dataset), each compared to a saved reference with a numeric tolerance. Implemented:
  `TestQSMMatrixPhantom.m` (methodQSMName), `TestBFRMatrixPhantom.m` (methodBFRName),
  `TestUnwrapMatrixPhantom.m` (methodUnwrapName). Each is parameterized over both
  method *and* matrix size (the same 5 sizes as `TestOddMatrixSize.m`) - only the
  baseline size (`[32 32 24]`) is compared numerically against a reference; the
  odd-size cases only check output shape/no-NaN, closing the gap `TestOddMatrixSize.m`
  documents (some pad/crop paths, e.g. Laplacian-based unwrap, only trigger for
  toolbox-dependent methods and are unreachable in Tier 1). Because the phantom is
  small, even iterative solvers finish quickly - this tier needs the proprietary
  toolboxes (MEDI/STI Suite/FANSI/SEGUE/MRITOOLS/HEIDI) installed locally, but **no
  real dataset**. Not run in CI. Run manually before tagging a release, or whenever
  you want fast toolbox-dependent coverage without a real dataset.
- **Tier 3 — real-dataset regression matrix** (`test/tier3_matrix/`): the same method
  matrix as Tier 2, but run once (no matrix-size parameterization) against a real
  dataset at its native, larger resolution - the slower, most realistic check.
  Implemented: `TestQSMMatrix.m`, `TestBFRMatrix.m`, `TestUnwrapMatrix.m`,
  `TestTwoPassMasking.m` (holds the QSM method fixed at FANSI and varies the
  two-pass mask refinement strategy instead — TKD/Direct Tikhonov apply the
  mask only after a closed-form k-space inversion, so two-pass masking has
  no effect on their output by construction; see
  sepia.documentation's Two-pass-masking.rst). Not yet implemented in any
  tier: R2*, a dedicated JSON-sidecar-contract audit. Not run in CI — it
  needs the proprietary toolboxes installed locally, and a real dataset.
  Run manually before tagging a release.

## Running Tier 1 locally

```
matlab -batch "cd('test/tools'); run_tier1"
```

No toolbox setup is required. `sepia_addpath` auto-creates a machine-local
`SpecifyToolboxesDirectory.m` (all paths empty) from the template on first call if
one doesn't already exist — this is expected and not an error.

A test row can show as:
- **Passed** — ran and matched its saved reference within tolerance.
- **Failed** — ran, but didn't match its reference (a real regression, or the
  reference is stale — see "Regenerating references" below).
- **Incomplete** (shown as "Filtered by assumption") — skipped because its
  reference file doesn't exist yet, or (Tier 2/3 only) its required toolbox/real
  dataset isn't available on this machine. This is *not* a failure.

`run_tier1.m` only fails the run (non-zero exit) on an actual **Failed** result.
`run_tier2.m`/`run_tier3.m` behave the same way.

## Running Tier 2 locally

On a machine with the relevant toolboxes installed and configured in
`SpecifyToolboxesDirectory.m` (see `getting_started/Installation.rst`), no real
dataset needed:

```
matlab -batch "cd('test/tools'); run_tier2"
```

Rows whose toolbox isn't installed will show as skipped (Incomplete), not failed —
that's expected on any machine that doesn't have every optional toolbox. Since the
phantom is small, this tier is fast even for iterative solvers - the full matrix
(method x matrix size) typically finishes in a few minutes.

## Running Tier 3 locally

Point `test/+sepiatest/get_real_dataset.m` at your data, either via two environment
variables:

```
export SEPIA_TEST_REAL_DATA_DIR=/path/to/your/dataset/input/directory
export SEPIA_TEST_REAL_DATA_MASK=/path/to/your/mask.nii.gz
```

or by writing them to `test/config/real_dataset.json` (gitignored, machine-local):

```json
{
  "inputDir": "/path/to/your/dataset/input/directory",
  "maskFile": "/path/to/your/mask.nii.gz"
}
```

`inputDir` is passed directly as `sepiaIO`'s `input` argument — either a directory
SEPIA can auto-detect (its own naming convention, or a BIDS-formatted directory with
multi-echo phase/magnitude NIfTI + JSON sidecars) — and `maskFile` as `maskFullName`;
the mask need not live inside `inputDir`. This has been validated against the
[QSM Consensus Paper example data](https://doi.org/10.1002/mrm.29048)'s
`converted/SIEMENS/Monopolar/GRE` directory (5-echo GRE, BIDS-style), using a brain
mask reused from that dataset's own previously-generated SEPIA derivatives.

Then, on a machine with the relevant toolboxes installed and configured in
`SpecifyToolboxesDirectory.m`:

```
matlab -batch "cd('test/tools'); run_tier3"
```

Rows whose toolbox isn't installed will show as skipped (Incomplete), not failed —
that's expected on any machine that doesn't have every optional toolbox. Tier 3 runs
on full-resolution real data, so expect each method to take on the order of a minute
to several minutes (closed-form methods like TKD are fast; iterative solvers like
MEDI/FANSI slower) — the full matrix can take tens of minutes. If you mainly want
toolbox-dependent method coverage without that wait, use Tier 2 instead.

## Interpreting a single failing row

Each `matlab.unittest` verification failure prints a clear message naming the exact
statistic and the actual vs. reference values, e.g.:

```
Tier1:QSM:TKD: stat "mean" outside tolerance (actual=-0.000456, reference=9.9995, relTol=0.0001, absTol=1e-06)
```

To re-run just one parameterized row, e.g. only the `iLSQR` smoke test:

```matlab
suite = matlab.unittest.TestSuite.fromClass(?TestSmokePhantom, 'ParameterName','qsmMethod','ParameterValue','iLSQR');
matlab.unittest.TestRunner.withTextOutput().run(suite)
```

## Regenerating references before a release

Do this only after an intentional change (an algorithmic fix, a MATLAB/toolbox
version bump, or a newly-added method needing its first baseline) — never as a
drive-by fix for a red test you don't understand.

```matlab
regenerate_reference('confirm', true, 'tier', 'tier1');                 % all Tier 1 references
regenerate_reference('confirm', true, 'tier', 'tier1', 'filter', 'TKD'); % just one

regenerate_reference('confirm', true, 'tier', 'tier2');                  % all Tier 2 references (synthetic phantom, all methods)
regenerate_reference('confirm', true, 'tier', 'tier2', 'filter', 'FANSI'); % just one (matches by method name across QSM/BFR/unwrap)

regenerate_reference('confirm', true, 'tier', 'tier3');                  % all Tier 3 references (real dataset - slow)
regenerate_reference('confirm', true, 'tier', 'tier3', 'filter', 'FANSI'); % just one
```

`regenerate_reference` refuses to run at all without `'confirm', true`, and prints an
old-vs-new percentage-change diff for every reference it touches so the change is
reviewable — treat a reference update as a reviewable change (`git diff` the
`.mat` file's `meta`/the printed diff, ideally with a second pair of eyes) before
committing it, exactly like any other code change. Tier 2/3 references are only
generated at the baseline matrix size (`[32 32 24]` for Tier 2's phantom) - the
odd-size cases in Tier 2 don't have (or need) their own reference.

## Known exclusions and caveats

- Phase unwrap **`'3D best path'`** is excluded entirely from automated comparison:
  it fills out-of-mask voxels with unseeded `rand()` and its own docstring says it
  only works on the DCCN cluster (it shells out to a compiled binary). Manual-only.
- **GPU paths** (NDI `isGPU=true`, QSMnet/LPCNN inference) are exercised CPU-only in
  this suite for determinism.
- **Tensor-MPPCA denoising** isn't exercised in the Tier 1 CI tier: its toolbox lives
  under `external/`, which is gitignored, so it won't exist on a fresh checkout (it
  auto-downloads on first use, or via `download_tMPPCA_toolbox.m`/`download_toolboxes.m`).
  Fine to include in a local Tier 2 or Tier 3 run.
- **Chi-separation / LPCNN / QSMnet+ / xQSM / BFRnet** (deep-learning methods) have no
  model/checkpoint files distributed with SEPIA or tracked in this repo, so this suite
  cannot exercise them at all yet — see `sepiatest.qsm_method_toolbox_key.m` /
  `sepiatest.bfr_method_toolbox_key.m` (the `'unavailable'` category).

## Adding a new method's regression row

Tier 2/3 test parameters are meant to be pulled from SEPIA's own method-name lists
(`methodUnwrapName`/`methodBFRName`/`methodQSMName`/... in
`configuration/sepia_configuration_*.m`, via `sepia_universal_variables`) rather than
hardcoded, so a newly-added method should automatically appear as a new
(initially-skipped, "no reference yet") row in both tiers. Add a tolerance-category
entry in `test/+sepiatest/tolerance_for_method.m` and run `regenerate_reference.m`
once per tier to create its baseline.
