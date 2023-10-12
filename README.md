# KBLTDARD

This is the data and code for the KBLTD.  Please cite if you use this code.

Data description:

"data_v2.mat" is the data format of matlab, which stores the HMDD v2.0 data set.

"data_v32.mat" is the data format of matlab, which stores the HMDD v3.2 data set.



Code description:

"KBLTD_ARD.m" code for the KBLTDARD method of this study.

"main_KBLTD_ARD_datav2.m" represents the code for performing 5-fold cross validation on HMDD v2.0.

"main_KBLTD_ARD_datav32.m" represents the code for performing 5-fold cross validation on HMDD v3.2.

"mianZ.m" is the main function. By running this function, you can get KBLTDARD's prediction results for HMDD v2.0 and HMDD v3.2 under CV_type and CV_triplet scenarios.

Before running "mianZ.m", download "tensor_toolbox-v3.1" and import it into the matlab running path.
