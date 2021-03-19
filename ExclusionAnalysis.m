function ExclusionAnalysis(Subregions,FolderPaths)
% Compares the results from WLS analyses with OLS analyses conducted using a range of exclusion threshold.
% INPUTS:
%   - Subregions: 'p1' and/or 'p2' depending on the tissue class of
% interest. Initialized in RunQUIQI 
%   - FolderPaths: structure containing the paths to the analysis
% folders - computed in PrepAnalysis
% Assumes the presence in each analysis folder of (only) ones T-statistical map of the mean effect across the cohort of interest. 
% 
% OUTPUTS: saved to disk. Maps of the difference in t-scores between WLS analyses and OLS analyses conducted 
%     - without excluding poor datasets ('Tchange_OLSToWLS.nii') 
%     - using the optimal exclusion threshold ('Tchange_ExclToWLS.nii').
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

Params=GetParams;
NMatlas=spm_read_vols(spm_vol(spm_select('FPList',Params.NMDir,'^label.*.nii$')));
FullCohortIndx=find(GetAnalParams('Exclusion').Exclusion==0);
OccurrenceIndx=zeros(size(Subregions,2),size(FolderPaths.CohortPaths,2),1);

TissueMask=ones(size(Subregions,2),size(NMatlas,1),size(NMatlas,2),size(NMatlas,3));
for tissuectr=1:size(Subregions,2)
    Mask=ones(size(NMatlas,1),size(NMatlas,2),size(NMatlas,3));
    Tmap=zeros(size(FolderPaths.CohortPaths,2),size(NMatlas,1),size(NMatlas,2),size(NMatlas,3));
    SavePath=fullfile(spm_str_manip(FolderPaths.CohortPaths{1},'h'),'AnalysisResults',Subregions{tissuectr});
    if ~exist(SavePath,'dir')
        mkdir(SavePath)
    end
    
    for repeatctr=1:size(FolderPaths.CohortPaths,2)
        CurrentPath=fullfile(FolderPaths.CohortPaths{repeatctr},FolderPaths.DataFolders{1},[Subregions{tissuectr} '_lambda_0']);
        TempMat=spm_read_vols(spm_vol(spm_select('FPList',CurrentPath,'^spmT.*.nii$')));
        Mask(find(TempMat==0))=0;%to ensure consistency of both explicit and implicit masks across exclusion thresholds
        
        Tmap(repeatctr,:,:,:)=spm_read_vols(spm_vol(spm_select('FPList',CurrentPath,'^spmT.*.nii$')));
    end
    
    CurrentPath=fullfile(FolderPaths.CohortPaths{1},FolderPaths.DataFolders{1},[Subregions{tissuectr} '_lambda_0']);
    Vsave=spm_vol(spm_select('FPList',CurrentPath,'^spmT.*.nii$'));
    
    [maxTmap,maxIndxmap]=max(Tmap,[],1);
    maxTmap=squeeze(maxTmap).*Mask;maxIndxmap=squeeze(maxIndxmap).*Mask;
    
    Vsave.fname=fullfile(SavePath,'maxT.nii');
    spm_write_vol(Vsave,maxTmap);
    Vsave.fname=fullfile(SavePath,'maxTIndx.nii');
    spm_write_vol(Vsave,maxIndxmap);
    for repeatctr=1:size(FolderPaths.CohortPaths,2)
        OccurrenceIndx(tissuectr,repeatctr)=size(find(maxIndxmap==repeatctr),1);
    end
    TissueMask(tissuectr,:,:,:)=Mask;
end
OccurrenceIndx=sum(OccurrenceIndx,1)/sum(sum(OccurrenceIndx,1),2)*100;
[~,maxOccurenceIndx]=max(sum(OccurrenceIndx,1),[],2);%Index of the exclusion fraction that yields the highest t-scores in OLS analysis, on average over all brain voxels.

save(fullfile(spm_str_manip(SavePath,'h'),'workspace.mat'), '-v7.3')

for tissuectr=1:size(Subregions,2)
    CurrentPath=fullfile(FolderPaths.CohortPaths{1},FolderPaths.DataFolders{1},[Subregions{tissuectr} '_lambda_0']);
    Vsave=spm_vol(spm_select('FPList',CurrentPath,'^spmT.*.nii$'));
    SavePath=fullfile(spm_str_manip(FolderPaths.CohortPaths{1},'h'),'AnalysisResults',Subregions{tissuectr});
    
    TFull_WLS=spm_read_vols(spm_vol(spm_select('FPList',fullfile(FolderPaths.CohortPaths{FullCohortIndx},'PDw_R2s',[Subregions{tissuectr} '_lambda_3']),'^spmT.*.nii$')));
    TFull_OLS=spm_read_vols(spm_vol(spm_select('FPList',fullfile(FolderPaths.CohortPaths{FullCohortIndx},'PDw_R2s',[Subregions{tissuectr} '_lambda_0']),'^spmT.*.nii$')));
    TExcl_OLS=spm_read_vols(spm_vol(spm_select('FPList',fullfile(FolderPaths.CohortPaths{maxOccurenceIndx},'PDw_R2s',[Subregions{tissuectr} '_lambda_0']),'^spmT.*.nii$')));
    
    Vsave.fname=fullfile(SavePath,'Tchange_ExclToWLS.nii');
    spm_write_vol(Vsave,(TFull_WLS-TExcl_OLS).*squeeze(TissueMask(tissuectr,:,:,:)));
    Vsave.fname=fullfile(SavePath,'Tchange_OLSToWLS.nii');
    spm_write_vol(Vsave,(TFull_WLS-TFull_OLS).*squeeze(TissueMask(tissuectr,:,:,:)));
    Vsave.fname=fullfile(SavePath,'AnalysisMask.nii');
    spm_write_vol(Vsave,squeeze(TissueMask(tissuectr,:,:,:)));    
end

end