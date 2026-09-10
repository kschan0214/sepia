# Two-pass masking exploration report

Fixed pipeline: ROMEO total field calculation (offset correction on) + VSHARP (STI suite) (radius=12, no polynomial refine) + FANSI (non-linear, TV, tol=0.1, maxiter=150, lambda=0.0001, mu1=0.01), isBET=1.

## SIEMENS Monopolar

### Strategy sanity check (lambda=0.7 / default threshold)

| strategy | volume ratio | dice | ran ok |
|---|---|---|---|
| Monoexponential decay model | 0.999 | 1.000 | 1 |
| Magnitude Gradient Field | 0.752 | 0.858 | 1 |
| Noise map | 0.968 | 0.984 | 1 |

### MGF lambda sweep

"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) to a voxel the refinement excluded - where a two-pass benefit should show up. "core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).

| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |
|---|---|---|---|---|---|---|
| 0.30 | 0.558 | 0.716 | 818183 | 0.0373 -> 0.0201 | [-0.0962,0.0919] -> [-0.0526,0.0509] | -0.0052 -> -0.0039 |
| 0.50 | 0.677 | 0.808 | 861919 | 0.0417 -> 0.0251 | [-0.1060,0.1076] -> [-0.0641,0.0663] | -0.0052 -> -0.0035 |
| 0.70 | 0.752 | 0.858 | 816174 | 0.0465 -> 0.0301 | [-0.1154,0.1257] -> [-0.0744,0.0832] | -0.0052 -> -0.0034 |
| 1.00 | 0.821 | 0.902 | 705362 | 0.0542 -> 0.0375 | [-0.1293,0.1572] -> [-0.0888,0.1101] | -0.0052 -> -0.0036 |
| 1.50 | 0.883 | 0.938 | 527620 | 0.0697 -> 0.0514 | [-0.1515,0.2339] -> [-0.1110,0.1626] | -0.0052 -> -0.0039 |
| 2.00 | 0.917 | 0.957 | 417434 | 0.0856 -> 0.0660 | [-0.1739,0.3098] -> [-0.1348,0.2263] | -0.0052 -> -0.0043 |
| 3.00 | 0.953 | 0.976 | 303256 | 0.1133 -> 0.0942 | [-0.2133,0.4128] -> [-0.1797,0.3419] | -0.0052 -> -0.0048 |

## SIEMENS Bipolar

### Strategy sanity check (lambda=0.7 / default threshold)

| strategy | volume ratio | dice | ran ok |
|---|---|---|---|
| Monoexponential decay model | 0.998 | 0.999 | 1 |
| Magnitude Gradient Field | 0.730 | 0.844 | 1 |
| Noise map | 0.962 | 0.981 | 1 |

### MGF lambda sweep

"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) to a voxel the refinement excluded - where a two-pass benefit should show up. "core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).

| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |
|---|---|---|---|---|---|---|
| 0.30 | 0.519 | 0.683 | 782506 | 0.0332 -> 0.0179 | [-0.0846,0.0842] -> [-0.0462,0.0437] | -0.0050 -> -0.0036 |
| 0.50 | 0.647 | 0.786 | 879735 | 0.0367 -> 0.0219 | [-0.0918,0.0978] -> [-0.0545,0.0572] | -0.0050 -> -0.0034 |
| 0.70 | 0.730 | 0.844 | 859121 | 0.0405 -> 0.0261 | [-0.0990,0.1128] -> [-0.0631,0.0722] | -0.0050 -> -0.0034 |
| 1.00 | 0.806 | 0.893 | 762532 | 0.0468 -> 0.0320 | [-0.1100,0.1392] -> [-0.0739,0.0953] | -0.0050 -> -0.0037 |
| 1.50 | 0.875 | 0.933 | 593694 | 0.0583 -> 0.0430 | [-0.1289,0.1902] -> [-0.0922,0.1395] | -0.0050 -> -0.0041 |
| 2.00 | 0.912 | 0.954 | 473678 | 0.0696 -> 0.0541 | [-0.1461,0.2451] -> [-0.1113,0.1839] | -0.0050 -> -0.0043 |
| 3.00 | 0.951 | 0.975 | 329759 | 0.0904 -> 0.0752 | [-0.1787,0.3323] -> [-0.1471,0.2730] | -0.0050 -> -0.0047 |

## GE Monopolar

### Strategy sanity check (lambda=0.7 / default threshold)

| strategy | volume ratio | dice | ran ok |
|---|---|---|---|
| Monoexponential decay model | 1.000 | 1.000 | 1 |
| Magnitude Gradient Field | 0.556 | 0.714 | 1 |
| Noise map | 0.987 | 0.994 | 1 |

### MGF lambda sweep

"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) to a voxel the refinement excluded - where a two-pass benefit should show up. "core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).

| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |
|---|---|---|---|---|---|---|
| 0.30 | 0.306 | 0.468 | 477454 | 0.0321 -> 0.0154 | [-0.0839,0.0787] -> [-0.0386,0.0394] | -0.0016 -> -0.0012 |
| 0.50 | 0.444 | 0.615 | 670816 | 0.0341 -> 0.0178 | [-0.0881,0.0843] -> [-0.0443,0.0444] | -0.0016 -> -0.0010 |
| 0.70 | 0.556 | 0.714 | 810539 | 0.0358 -> 0.0202 | [-0.0920,0.0899] -> [-0.0497,0.0509] | -0.0016 -> -0.0008 |
| 1.00 | 0.676 | 0.806 | 883893 | 0.0390 -> 0.0243 | [-0.0992,0.1007] -> [-0.0595,0.0639] | -0.0016 -> -0.0007 |
| 1.50 | 0.792 | 0.884 | 822456 | 0.0450 -> 0.0310 | [-0.1113,0.1241] -> [-0.0740,0.0867] | -0.0016 -> -0.0010 |
| 2.00 | 0.856 | 0.922 | 700936 | 0.0514 -> 0.0379 | [-0.1245,0.1494] -> [-0.0887,0.1128] | -0.0016 -> -0.0012 |
| 3.00 | 0.919 | 0.958 | 494148 | 0.0648 -> 0.0515 | [-0.1514,0.2049] -> [-0.1159,0.1635] | -0.0016 -> -0.0014 |

## GE Bipolar

### Strategy sanity check (lambda=0.7 / default threshold)

| strategy | volume ratio | dice | ran ok |
|---|---|---|---|
| Monoexponential decay model | 0.999 | 0.999 | 1 |
| Magnitude Gradient Field | 0.596 | 0.747 | 1 |
| Noise map | 0.987 | 0.994 | 1 |

### MGF lambda sweep

"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) to a voxel the refinement excluded - where a two-pass benefit should show up. "core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).

| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |
|---|---|---|---|---|---|---|
| 0.30 | 0.343 | 0.511 | 525770 | 0.0354 -> 0.0163 | [-0.0933,0.0882] -> [-0.0420,0.0416] | -0.0023 -> -0.0018 |
| 0.50 | 0.486 | 0.654 | 723785 | 0.0378 -> 0.0191 | [-0.0988,0.0944] -> [-0.0493,0.0473] | -0.0023 -> -0.0016 |
| 0.70 | 0.596 | 0.747 | 847258 | 0.0399 -> 0.0218 | [-0.1046,0.1002] -> [-0.0564,0.0547] | -0.0023 -> -0.0014 |
| 1.00 | 0.708 | 0.829 | 880162 | 0.0437 -> 0.0265 | [-0.1131,0.1119] -> [-0.0672,0.0680] | -0.0023 -> -0.0015 |
| 1.50 | 0.814 | 0.897 | 789026 | 0.0504 -> 0.0339 | [-0.1283,0.1354] -> [-0.0846,0.0927] | -0.0023 -> -0.0016 |
| 2.00 | 0.871 | 0.931 | 658304 | 0.0577 -> 0.0415 | [-0.1442,0.1626] -> [-0.1015,0.1195] | -0.0023 -> -0.0017 |
| 3.00 | 0.928 | 0.963 | 452550 | 0.0728 -> 0.0568 | [-0.1748,0.2288] -> [-0.1342,0.1756] | -0.0023 -> -0.0020 |

## PHILIPS Monopolar CLEAR

### Strategy sanity check (lambda=0.7 / default threshold)

| strategy | volume ratio | dice | ran ok |
|---|---|---|---|
| Monoexponential decay model | 1.000 | 1.000 | 1 |
| Magnitude Gradient Field | 0.589 | 0.741 | 1 |
| Noise map | 0.994 | 0.997 | 1 |

### MGF lambda sweep

"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) to a voxel the refinement excluded - where a two-pass benefit should show up. "core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).

| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |
|---|---|---|---|---|---|---|
| 0.30 | 0.349 | 0.517 | 503271 | 0.0371 -> 0.0176 | [-0.0965,0.0917] -> [-0.0471,0.0434] | -0.0044 -> -0.0034 |
| 0.50 | 0.487 | 0.655 | 674526 | 0.0397 -> 0.0209 | [-0.1032,0.0990] -> [-0.0554,0.0524] | -0.0044 -> -0.0029 |
| 0.70 | 0.589 | 0.741 | 768821 | 0.0424 -> 0.0243 | [-0.1098,0.1063] -> [-0.0647,0.0619] | -0.0044 -> -0.0024 |
| 1.00 | 0.693 | 0.819 | 799967 | 0.0467 -> 0.0293 | [-0.1206,0.1204] -> [-0.0773,0.0765] | -0.0044 -> -0.0020 |
| 1.50 | 0.795 | 0.886 | 732515 | 0.0544 -> 0.0377 | [-0.1386,0.1477] -> [-0.0977,0.1048] | -0.0044 -> -0.0019 |
| 2.00 | 0.854 | 0.921 | 623509 | 0.0632 -> 0.0465 | [-0.1565,0.1812] -> [-0.1182,0.1346] | -0.0044 -> -0.0021 |
| 3.00 | 0.913 | 0.955 | 450685 | 0.0804 -> 0.0619 | [-0.1896,0.2584] -> [-0.1486,0.1914] | -0.0044 -> -0.0026 |

## PHILIPS Bipolar CLEAR

### Strategy sanity check (lambda=0.7 / default threshold)

| strategy | volume ratio | dice | ran ok |
|---|---|---|---|
| Monoexponential decay model | 0.999 | 0.999 | 1 |
| Magnitude Gradient Field | 0.618 | 0.764 | 1 |
| Noise map | 0.994 | 0.997 | 1 |

### MGF lambda sweep

"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) to a voxel the refinement excluded - where a two-pass benefit should show up. "core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).

| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |
|---|---|---|---|---|---|---|
| 0.30 | 0.370 | 0.540 | 531469 | 0.0360 -> 0.0179 | [-0.0939,0.0864] -> [-0.0473,0.0437] | -0.0055 -> -0.0040 |
| 0.50 | 0.513 | 0.679 | 710371 | 0.0386 -> 0.0211 | [-0.1002,0.0948] -> [-0.0549,0.0518] | -0.0055 -> -0.0035 |
| 0.70 | 0.618 | 0.764 | 795246 | 0.0414 -> 0.0245 | [-0.1068,0.1045] -> [-0.0629,0.0623] | -0.0055 -> -0.0030 |
| 1.00 | 0.720 | 0.837 | 797928 | 0.0459 -> 0.0293 | [-0.1173,0.1208] -> [-0.0740,0.0789] | -0.0055 -> -0.0026 |
| 1.50 | 0.816 | 0.898 | 695190 | 0.0542 -> 0.0372 | [-0.1345,0.1529] -> [-0.0903,0.1099] | -0.0055 -> -0.0027 |
| 2.00 | 0.868 | 0.929 | 572988 | 0.0635 -> 0.0463 | [-0.1528,0.1902] -> [-0.1095,0.1435] | -0.0055 -> -0.0029 |
| 3.00 | 0.921 | 0.959 | 409275 | 0.0818 -> 0.0638 | [-0.1875,0.2759] -> [-0.1441,0.2092] | -0.0055 -> -0.0036 |

## PHILIPS Monopolar SYNERGY

### Strategy sanity check (lambda=0.7 / default threshold)

| strategy | volume ratio | dice | ran ok |
|---|---|---|---|
| Monoexponential decay model | 1.000 | 1.000 | 1 |
| Magnitude Gradient Field | 0.577 | 0.732 | 1 |
| Noise map | 0.996 | 0.998 | 1 |

### MGF lambda sweep

"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) to a voxel the refinement excluded - where a two-pass benefit should show up. "core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).

| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |
|---|---|---|---|---|---|---|
| 0.30 | 0.326 | 0.492 | 466476 | 0.0328 -> 0.0169 | [-0.0852,0.0801] -> [-0.0442,0.0419] | -0.0042 -> -0.0033 |
| 0.50 | 0.466 | 0.636 | 647402 | 0.0347 -> 0.0197 | [-0.0906,0.0859] -> [-0.0508,0.0484] | -0.0042 -> -0.0029 |
| 0.70 | 0.577 | 0.732 | 767809 | 0.0366 -> 0.0225 | [-0.0949,0.0923] -> [-0.0577,0.0562] | -0.0042 -> -0.0025 |
| 1.00 | 0.691 | 0.817 | 813702 | 0.0398 -> 0.0259 | [-0.1022,0.1041] -> [-0.0661,0.0676] | -0.0042 -> -0.0022 |
| 1.50 | 0.798 | 0.887 | 743998 | 0.0460 -> 0.0307 | [-0.1146,0.1279] -> [-0.0767,0.0887] | -0.0042 -> -0.0020 |
| 2.00 | 0.856 | 0.923 | 633294 | 0.0526 -> 0.0376 | [-0.1264,0.1536] -> [-0.0910,0.1138] | -0.0042 -> -0.0022 |
| 3.00 | 0.915 | 0.956 | 445209 | 0.0675 -> 0.0520 | [-0.1506,0.2242] -> [-0.1165,0.1678] | -0.0042 -> -0.0027 |

## PHILIPS Bipolar SYNERGY

### Strategy sanity check (lambda=0.7 / default threshold)

| strategy | volume ratio | dice | ran ok |
|---|---|---|---|
| Monoexponential decay model | 0.999 | 0.999 | 1 |
| Magnitude Gradient Field | 0.595 | 0.746 | 1 |
| Noise map | 0.995 | 0.998 | 1 |

### MGF lambda sweep

"boundary" = voxels retained in the refined mask but adjacent (within 3 vox) to a voxel the refinement excluded - where a two-pass benefit should show up. "core" = mask eroded by 9 vox, a deep-brain stability check (should stay ~unchanged).

| lambda | mask vol ratio | dice | n boundary vox | boundary std (pass1->combined) | boundary p1/p99 (pass1->combined) | core mean (pass1->combined) |
|---|---|---|---|---|---|---|
| 0.30 | 0.342 | 0.510 | 487392 | 0.0339 -> 0.0168 | [-0.0881,0.0856] -> [-0.0448,0.0417] | -0.0034 -> -0.0027 |
| 0.50 | 0.486 | 0.654 | 670254 | 0.0362 -> 0.0199 | [-0.0943,0.0919] -> [-0.0525,0.0490] | -0.0034 -> -0.0024 |
| 0.70 | 0.595 | 0.746 | 776701 | 0.0384 -> 0.0230 | [-0.1004,0.0980] -> [-0.0597,0.0582] | -0.0034 -> -0.0021 |
| 1.00 | 0.705 | 0.827 | 810040 | 0.0420 -> 0.0276 | [-0.1093,0.1096] -> [-0.0712,0.0732] | -0.0034 -> -0.0017 |
| 1.50 | 0.808 | 0.894 | 731649 | 0.0488 -> 0.0348 | [-0.1243,0.1329] -> [-0.0889,0.0994] | -0.0034 -> -0.0017 |
| 2.00 | 0.864 | 0.927 | 613772 | 0.0560 -> 0.0414 | [-0.1388,0.1593] -> [-0.1032,0.1229] | -0.0034 -> -0.0018 |
| 3.00 | 0.921 | 0.959 | 424583 | 0.0711 -> 0.0559 | [-0.1652,0.2248] -> [-0.1324,0.1755] | -0.0034 -> -0.0022 |

