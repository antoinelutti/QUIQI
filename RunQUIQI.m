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
% 2. 'Inserting an index of image quality into the analysis of MRI data',
% Lutti et al, in prep.

Params=GetParams;
eval(['load ' fullfile(Params.HomeDir,'Subject_Details.mat')]);

lambda={[0 1 2 3],0,1,2,3,4,5};
Subregions={'p1','p2'};DataType={'PDw_R2s','MTw_R2s','T1w_R2s'};AnalType='Full';
% Includes GM regional analysis:
% lambda={0,1,2,3,4,5};
% Subregions=cat(2,{'p1','p2'},getROIpairs('GM'));DataType={'PDw_R2s'};AnalType='Full';
[QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
RunAnalysis(QUIQI,AnalType);
FreeEnergyAnalysis(lambda,Subregions,FolderPaths)

lambda={0,3};
Subregions={'p1','p2'};DataType={'PDw_R2s'};AnalType='Residuals';
[QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
RunAnalysis(QUIQI,AnalType);
MDIvsResAnalysis(QUIQI,FolderPaths);

lambda={0,3};%The 'Exclusion'analysis is also conducted for lambda~=0 for simplicity, although the results are only analysed for the size of the excluded datasets is zero.
Subregions=cat(2,{'p1','p2'});DataType={'PDw_R2s'};AnalType='Exclusion';
[QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
RunAnalysis(QUIQI,AnalType);
ExclusionAnalysis(Subregions,FolderPaths)

lambda={0,3};
Subregions={'p1','p2'};DataType={'PDw_R2s'};AnalType='Specificity';
[QUIQI,FolderPaths]=PrepAnalysis(Subject_Details,lambda,DataType,Subregions,AnalType);
RunAnalysis(QUIQI,AnalType);
FalsePositiveCount(FolderPaths,'^spmF.*.(img|nii)$',[15.07 15.28 13.09 13.19],[7.33 7.33 7.33 7.33]);


end



