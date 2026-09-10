# Chi-separation QSM Addon

This addon provides χ-separation (Chi-separation) QSM algorithms (e.g. QSMnet+, xsepnet, R2PRIME-net) for SEPIA, using ONNX network checkpoints.

## 1. Download the checkpoint files

The official source for the χ-separation toolbox is the SNU-LIST repository: https://github.com/SNU-LIST/chi-separation

There is no direct download link or GitHub release - the repository requires submitting a Google Form, after which a download link for the toolbox is sent by email:

**https://forms.gle/nhJahF86zpMgEKvM9**

Unpacking the toolbox gives a `Chisep_Toolbox_v1.2/` folder containing (among other things) a `models/` subfolder. Copy the ONNX model checkpoint files (`240531_R2PRIMEnet.onnx`, `240904_QSMnet.onnx`, `240904_xsepnet.onnx`, `R2PNET_7T.onnx`) and their associated normalization factor files (`norm_factor.mat`, `xsepnet_train_patch_norm_factor_inplane_largedegree_romeo_arlo.mat`) from that `models/` subfolder into this addon's own [models/](models/) folder.

## 2. Set `ChiSepNet_HOME`

`ChiSepNet_HOME` is configured centrally in `SpecifyToolboxesDirectory.m` (like every other optional toolbox), rather than hand-edited in this addon's own files. Set it to the path of this addon folder (i.e. `addons/qsm/Chi-separation`, where the `models/` folder from step 1 lives), either by editing `SpecifyToolboxesDirectory.m` directly:

```matlab
ChiSepNet_HOME = '/path/to/your/sepia/addons/qsm/Chi-separation';
```

or via the GUI's Utility tab → Manage Dependency panel.

## 3. Required MATLAB Add-On

Loading the ONNX models requires `importONNXNetwork`, which needs the **Deep Learning Toolbox Converter for ONNX Model Format** support package.

To install it, open MATLAB's **Add-On Explorer** and search for "Deep Learning Toolbox Converter for ONNX Model Format", then install it.
