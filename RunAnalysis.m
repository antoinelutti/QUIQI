function RunAnalysis(QUIQI,AnalType)
% Runs image analysis.
% INPUTS: 
%     - QUIQI: structure containing the information required for the analysis. Computed by PrepAnalysis.m
%     - AnalType: analysis type. Initialized in RunQUIQI.m.
%
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

for subsetctr=1:size(QUIQI,2)% Analysis folder is deleted if already exists
    if exist(fullfile(char(QUIQI(subsetctr).CohortPath),char(QUIQI(subsetctr).AnalDir),'SPM.mat'),'file')
        s=['delete ' fullfile(char(QUIQI(subsetctr).CohortPath),char(QUIQI(subsetctr).AnalDir),'SPM.mat')];eval(s);
    end
    if exist(fullfile(char(QUIQI(subsetctr).CohortPath),char(QUIQI(subsetctr).AnalDir)))==7
        rmdir(fullfile(char(QUIQI(subsetctr).CohortPath),char(QUIQI(subsetctr).AnalDir)),'s');
    end
end
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool(12);
parfor subsetctr=1:size(QUIQI,2)
    Subject_Details_subset=load(fullfile(char(QUIQI(subsetctr).CohortPath),'Subject_Details.mat'));
    Subject_Details_subset=Subject_Details_subset.Subject_Details;
    
    CurrentDir=fullfile(char(QUIQI(subsetctr).CohortPath),char(QUIQI(subsetctr).AnalDir));
    if(~exist(CurrentDir,'dir'))
        mkdir(CurrentDir)
    end
    copyfile(spm_select('FPList',char(QUIQI(subsetctr).CohortPath),['^ExplicitMask_' QUIQI(subsetctr).ROI '.nii']),CurrentDir);% copy masks to local folder
    myparsave(fullfile(CurrentDir,'ReMLcomps'),QUIQI(subsetctr).ReML,'ReMLcomps');
    
    Covariates=[Subject_Details_subset.confound];
    DM=MakeDM(Subject_Details_subset,QUIQI(subsetctr).MotionReg{1},AnalType);
    RunSPMAnal(Subject_Details_subset,CurrentDir,QUIQI(subsetctr).InputData,QUIQI(subsetctr).ROI,QUIQI(subsetctr).ReML,[],'explicit',AnalType,DM,[Covariates.BrainVol]')
end

end

function RunSPMAnal(Subject_Details,CurrentDir,InputData,ROI,CoVReml,Weights,masking,AnalType,DM,BrainVol)

Params=GetParams;
AnalParams=GetAnalParams(AnalType);
ExplicitMask=spm_select('FPList',CurrentDir,['^ExplicitMask_' ROI '.nii']);

if ~isempty(Weights)
    save(fullfile(CurrentDir,'Weights'),'Weights')
end

fileID = fopen(fullfile(CurrentDir,'Scans.txt'),'w');
for ctr=1:size(Subject_Details,2)
    P=spm_select('FPList',fullfile(Params.DataDir,char(Subject_Details(ctr).ID),Params.DataSubDir),['^' InputData]);
    fprintf(fileID,'%s\n',spm_str_manip(P,'t'));
    Scans{ctr}=P;
end
fclose(fileID);

Scans=Scans';

spm_jobman('initcfg');
clear matlabbatch
matlabbatch{1}.spm.stats.factorial_design.dir = {CurrentDir};
%%
if strcmp(AnalType,'GroupComparison')
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = Scans(1:AnalParams.NGroup1);
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = Scans(AnalParams.NGroup1+1:end);
    matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;    
else
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = Scans;
end
%%
%%
for ctr=1:size(DM.mat,2)
    matlabbatch{1}.spm.stats.factorial_design.cov(ctr).c = DM.mat(:,ctr);
    %%
    matlabbatch{1}.spm.stats.factorial_design.cov(ctr).cname = DM.text{ctr};
    matlabbatch{1}.spm.stats.factorial_design.cov(ctr).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(ctr).iCC = 1;
end

matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
if strcmp('explicit',masking)
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {ExplicitMask};
elseif strcmp('implicit',masking)
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
elseif strcmp('bothplicit',masking)
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {ExplicitMask};
end
%%
matlabbatch{1}.spm.stats.factorial_design.globalc.g_user.global_uval = BrainVol';
%%
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 3;
spm_jobman('run', matlabbatch);
   
if ~isempty(CoVReml)
    load(fullfile(CurrentDir,'SPM.mat'))
    if isfield(SPM.xVi,'Vi')% Group comparison
        GroupIndx={};DataSize=size(SPM.xVi.Vi{1},1);
        for ctr=1:size(SPM.xVi.Vi,2)
            GroupIndx{ctr}=find(diag(SPM.xVi.Vi{ctr})~=0);
        end
    else
        DataSize=size(SPM.xVi.V,1);
        GroupIndx{1}=find(diag(SPM.xVi.V)~=0);
    end
    SPM=rmfield(SPM,'xVi');
    ctr1d=0;
    for ctr=1:size(GroupIndx,2)
        for ctr2=1:size(CoVReml,2)
            ctr1d=ctr1d+1;
            DiagTerms=zeros(DataSize,1);DiagTerms(GroupIndx{ctr},1)=CoVReml(GroupIndx{ctr},ctr2);
            SPM.xVi.Vi(ctr1d)={sparse(diag(DiagTerms))};
        end
    end
    if isfield(DM,'ReMLFcontrast')
        SPM.xVi.Fcontrast = DM.ReMLFcontrast;
    end
    save(fullfile(CurrentDir,'SPM.mat'), 'SPM')
end
if ~isempty(Weights)&&~isequal(Weights,ones(size(Weights)))
    load(fullfile(CurrentDir,'SPM.mat'))
    SPM.xX.K=sparse(diag(Weights));
    save(fullfile(CurrentDir,'SPM.mat'), 'SPM')
end

spm_jobman('initcfg');
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(CurrentDir,'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = AnalParams.SaveResiduals;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch);

if strcmp(AnalType,'GroupComparison')
    spm_jobman('initcfg');
    clear matlabbatch
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(CurrentDir,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Tests > Controls';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Controls > Tests';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run', matlabbatch);
% elseif strcmp(AnalType,'Exclusion')
%     spm_jobman('initcfg');
%     clear matlabbatch
%     matlabbatch{1}.spm.stats.con.spmmat = {fullfile(CurrentDir,'SPM.mat')};
%     matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Mean';
%     matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 0];
%     matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
%     matlabbatch{1}.spm.stats.con.delete = 1;
%     spm_jobman('run', matlabbatch);
else
    spm_jobman('initcfg');
    clear matlabbatch
    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(CurrentDir,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = DM.desc;
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = cat(2,zeros(DM.size,1),eye(DM.size));
    matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    spm_jobman('run', matlabbatch);
end

if strcmp(AnalType,'GroupComparison')||strcmp(AnalType,'Specificity')
    ToBedeleted=spm_select('FPList',CurrentDir,'^(ExplicitMask|ReMLMask|beta|ResMS|RPV|con|mask|ess).*.(mat|img|nii)$');
    for ctr=1:size(ToBedeleted,1)
        delete(deblank(ToBedeleted(ctr,:)));
    end
end

end


function myparsave(SaveDir, ToBeSaved,SaveName)
eval([SaveName '= ToBeSaved;'])
save(SaveDir,SaveName, '-v7.3')
end