function RunQUIQI
% Main function enabling the insertion of a motion degradation index into
% the analysis of MRI data. 
%
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland
% 
% REFERENCES
% 1. 'QUIQI – using a QUality Index for the analysis of Quantitative Imaging
% data'. Di Domenicantonio et al., Proceedings of the Annual Meeting of the
% International Society for Magnetic Resonance in Medicine, 2020.
%  2. 'Inserting an index of image quality into the statistical analysis of MRI data',
% Lutti et al, in prep.


Params=GetParams;
eval(['load ' fullfile(Params.DataDir,'Subject_Details.mat')]);

lambda={0,[0 1 2 3],0,1,2,3,4,5};
Subregions={'p1','p2'};DataType={'PDw_R2s'};AnalType='Full';
% lambda={0,1,2,3,4,5};% Commented due to size of output data (regional analysis)
% Subregions=cat(2,{'p1','p2'},getROIpairs('GM'));DataType={'PDw_R2s'};AnalType='Full';
[QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
RunAnalysis(QUIQI,AnalType);
% !! The following requires the versions of spm_spm.m and spm_est_non_shericity.m
% provided in this package!!
FreeEnergyAnalysis(lambda,Subregions,FolderPaths)

lambda={0,3};
Subregions={'p1','p2'};DataType={'PDw_R2s'};AnalType='Residuals';
[QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
RunAnalysis(QUIQI,AnalType);
% !! The following requires the versions of spm_spm.m and spm_est_non_shericity.m
% provided in this package!!
MDIvsResAnalysis(QUIQI,FolderPaths);

lambda={0};
Subregions=cat(2,{'p1','p2'});DataType={'PDw_R2s'};AnalType='Exclusion';
[QUIQI,~]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
RunAnalysis(QUIQI,AnalType);

% lambda={0,3};% Commented because requires data within a narrow age-range (see README.pdf)
% Subregions={'p1','p2'};DataType={'PDw_R2s'};AnalType='MotionBias';
% [QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
% RunAnalysis(QUIQI,AnalType);
% MDIvsResAnalysis(QUIQI,FolderPaths);
% 
% lambda={0,3};% Commented because requires data within a narrow age-range (see README.pdf)
% Subregions={'p1','p2'};DataType={'PDw_R2s'};AnalType='GroupComparison';
% [QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
% RunAnalysis(QUIQI,AnalType);
% FalsePositiveCount(FolderPaths,'^spmT.*.(img|nii)$',[4.83 4.89 4.47 4.53],[3.16 3.16 3.16 3.16]);

lambda={0,3};
Subregions={'p1','p2'};DataType={'PDw_R2s'};AnalType='Specificity';
[QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
RunAnalysis(QUIQI,AnalType);
FalsePositiveCount(FolderPaths,'^spmF.*.(img|nii)$',[15.07 15.28 13.09 13.19],[7.33 7.33 7.33 7.33]);

end



