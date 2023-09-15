%%%%%%%%%%%%%
clc
clear

% % %%%%%%%%%%%%%%%%%%%%%%%%HMDD v2.0
load('data_v2.mat')
main_KBLTD_ARD_datav2(data_v2,'KBLTD_ARD','cv_type');    %%%%CV_type
main_KBLTD_ARD_datav2(data_v2,'KBLTD_ARD','cv_triplet');   %%%%CV_triplet



% % %%%%%%%%%%%%%%%%%%%%%%%%HMDD v3.2
load('data_v32.mat')
main_KBLTD_ARD_datav32(data_v32,'KBLTD_ARD','cv_type');    %%%%CV_type
main_KBLTD_ARD_datav32(data_v32,'KBLTD_ARD','cv_triplet');   %%%%CV_triplet



