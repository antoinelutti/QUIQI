function PrepareMasksAndCohort(CurrentPath,CreateNewMask,Subregions)
% Creates explicit mask for image analysis.
% INPUTS:
%     - CurrentPath: location of the analysis folder
%     - CreateNewMask: binary flag. If CreateNewMask is set to 1, a new explicit mask will be computed
%     even if one already exists.
%     - Subregions: region of interest included in the computed explicit
%     mask. Initialized in RunQUIQI.m.

% This function creates an explicit mask if it's non-existent at the
% location specified by CurrentPath or if CreateNewMask is set to 1 (e.g.
% because a new demographic information file was created at the
% start of the analysis, suggesting a new analysis).
% The mask will be defined over a region defined by Subregions. Values for Subregions might be p1/p2
% for grey/white matter or cortical regions labels of the neuromorphometric atlas.
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

MakeNewMask=0;
for ctr=1:size(Subregions,2)
    RegionStr=RegionLabeltoName(Subregions(ctr));
    if (~exist(fullfile(CurrentPath,['ExplicitMask_' char(RegionStr{1}) '.nii']),'file') || CreateNewMask==1)
        MakeNewMask=1;
    end
end
if MakeNewMask==1
    load(fullfile(CurrentPath,'Subject_Details.mat'));
    MakeExplicitMask(CurrentPath,Subject_Details,Subregions,RegionLabeltoName(Subregions));
end
end

function MakeExplicitMask(SavePath,Subject_Details,Subregions,RegionStr)

Params=GetParams;
c1str='^mwc1.*.(img|nii)$';
c2str='^mwc2.*.(img|nii)$';

for ctrsubj=1:size(Subject_Details,2)
    segfolder=char(fullfile(Params.DataDir,char(Subject_Details(ctrsubj).ID),Params.DataSubDir));
    c1=spm_read_vols(spm_vol(spm_select('FPList',segfolder,c1str)));
    c2=spm_read_vols(spm_vol(spm_select('FPList',segfolder,c2str)));
    if ctrsubj==1
        Vsave=spm_vol(spm_select('FPList',segfolder,c1str));
        sc1=single(zeros(Vsave.dim(1),Vsave.dim(2),Vsave.dim(3)));
        sc2=single(zeros(Vsave.dim(1),Vsave.dim(2),Vsave.dim(3)));
    end
    sc1=sc1+single(smooth3(c1,'gaussian',[3 3 3]));
    sc2=sc2+single(smooth3(c2,'gaussian',[3 3 3]));
end
sc1=sc1/size(Subject_Details,2);sc2=sc2/size(Subject_Details,2);
GMmask=zeros(size(sc1));WMmask=zeros(size(sc2));
GMmask(find(sc1>sc2&(sc1+sc2)>0.5))=1;WMmask(find(sc2>sc1&(sc1+sc2)>0.5))=1;
if ~exist(SavePath,'dir')
    mkdir(SavePath)
end
Vsave.fname=fullfile(SavePath,'ExplicitMask_p1.nii');spm_write_vol(Vsave,GMmask);
Vsave.fname=fullfile(SavePath,'ExplicitMask_p2.nii');spm_write_vol(Vsave,WMmask);

NMatlas=spm_read_vols(spm_vol(spm_select('FPList',Params.NMDir,'^label.*.nii$')));

for ctr=1:size(Subregions,2)
    if ~strcmp(Subregions{ctr},'p1') & ~strcmp(Subregions{ctr},'p2')
        LocalMask=zeros(size(GMmask));
        if ismember(Subregions{ctr},Params.BrainRegions.GMregions)==1
            LocalMask(ismember(NMatlas,Subregions{ctr}) & GMmask == 1)=1;
        elseif ismember(Subregions{ctr},Params.BrainRegions.WMregions)==1
            LocalMask(ismember(NMatlas,Subregions{ctr}) & WMmask == 1)=1;
        end
        Vsave.fname=fullfile(SavePath,['ExplicitMask_' RegionStr{ctr} '.nii']);
        spm_write_vol(Vsave,LocalMask);        
    end    
end
end

