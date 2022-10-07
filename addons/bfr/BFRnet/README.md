# Set-up BFRnet in SEPIA

1. Download [deepMRI](https://github.com/sunhongfu/deepMRI) from GitHub
2. Download the pre-trained BFRnet using this [link](https://www.dropbox.com/sh/q678oapc65evrfa/AADh2CGeUzhHh6q9t3Fe3fVVa?dl=0) (see https://github.com/sunhongfu/deepMRI/tree/master/BFRnet)
3. Specify the full path to deepMRI code as 'deepMRI_HOME' in setup_BFRnet_environment.m 
4. Specify the full path to the pre-trained network (should be checkpoints/BFRnet_L2_64PS_24BS_45Epo_NewHCmix.mat) from (2) as 'checkpoints' in setup_BFRnet_environment.m 
