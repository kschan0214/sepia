# Chi-separation QSM Addon

This addon provides χ-separation (Chi-separation) QSM algorithms (e.g. QSMnet+, xsepnet, R2PRIME-net) for SEPIA, using ONNX network checkpoints.

## 1. Download the checkpoint files

The ONNX model checkpoint files (e.g. `240531_R2PRIMEnet.onnx`, `240904_QSMnet.onnx`, `240904_xsepnet.onnx`, `R2PNET_7T.onnx`) and their associated normalization factor files (`norm_factor.mat`, `xsepnet_train_patch_norm_factor_inplane_largedegree_romeo_arlo.mat`) need to be placed in the [models/](models/) folder.

Download link: **TODO — to be provided**

## 2. Update the path in `setup_Chi_sepnet_environment.m`

Open [setup_Chi_sepnet_environment.m](setup_Chi_sepnet_environment.m) and update `home_directory` to point to the actual location of this addon (i.e. the `Chi-separation` folder) on your system, e.g.:

```matlab
home_directory = '/path/to/your/sepia/addons/qsm/Chi-separation';
```

## 3. Required MATLAB Add-On

Loading the ONNX models requires `importONNXNetwork`, which needs the **Deep Learning Toolbox Converter for ONNX Model Format** support package.

To install it, open MATLAB's **Add-On Explorer** and search for "Deep Learning Toolbox Converter for ONNX Model Format", then install it.
