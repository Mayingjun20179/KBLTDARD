Matlab code of Kernel Bayesian logistic tensor decomposition with automatic rank determination for predicting multiple types of miRNA-disease associations
(Written by Yingjun Ma 2024)


To run the code:
1. Change Matlab work directory to "/KBLTDARD/".
2. Run  "loadpah" code to add the current folder and subfolders into Matlab path searching list.
3. Open and run the demo file. 


Demo code:   mianZ.m

The “Cv_experiment” folder contains related experiments on both datasets
I.  main_KBLTD_ARD_datav2.m:    Performing 5-fold cross validation on HMDD v2.0
II. main_KBLTD_ARD_datav32.m   Performing 5-fold cross validation on HMDD v3.2


The “Algorithms” folder contains the relevant calculation code for KBLTD_ARD:
KBLTD_ARD.m：Kernel Bayesian logistic tensor decomposition with automatic rank determination.

The Dataset folder contains all experimental data in this paper:
I. data_v2.mat：Benchmark dataset for HMDD v2.0
II.data_v32.mat： Benchmark dataset for HMDD v3.2



In this package, we used the tensor tensor_toolbox-v3.1, which is downloaded from (https://gitlab.com/tensors/tensor_toolbox/-/releases/v3.1)
