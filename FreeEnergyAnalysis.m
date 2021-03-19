function FreeEnergyAnalysis(lambda,Subregions,FolderPaths)
% Analysis of the Free Energy change estimated from ReML for different noise models
% INPUTS: 
%     - lambda: cell array of power values of the MDI, initialized in RunQUIQI.m
%     - Subregions:cell array of regions  of interests for the analysis. Initialized in RunQUIQI.m.
%     - FolderPaths: structure of analysis folder pathts. Computed by PrepAnalysis.m
% OUTPUTS: saved to disk in 'ReMLAnal' folder
%     - For analysis over a whole tissue class (Subregions = p1 or p2), the output of the analysis is a matlab figure containing a histogram of free energy values.
%     - For local analysis (when Subregions contains atlas regional labels), the outputs are maps of Free Energy difference at the regional level. 
%       These differences are calculated between OLS analyses, WLS analyses yielding the maximum local Free Energy and WLS analyses
%       from a reference noise model defined by Params.ReMLAnal.RefLambda.
%
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

Params=GetParams;
RootPath=FolderPaths.CohortPaths{1};
LambdaRefIndx=1;OLSIndx=1;ManualExclusionIndx=[];
for lambdactr=1:size(lambda,2)
    if isnumeric(lambda{lambdactr})
        if isequal(lambda{lambdactr},Params.ReMLAnal.RefLambda)
            LambdaRefIndx=lambdactr;
        elseif lambda{lambdactr}==0
            OLSIndx=lambdactr;
        end
    else
        if isempty(ManualExclusionIndx)
            ManualExclusionIndx=lambdactr;
        else
            ManualExclusionIndx=cat(2,ManualExclusionIndx,lambdactr);
        end
    end
end
NotOLS=find(ismember(linspace(1,size(lambda,2),size(lambda,2)),OLSIndx)==0);


GMIndx=find(strcmp(repmat({'p1'},[1,size(Subregions,2)]),Subregions)==1);
WMIndx=find(strcmp(repmat({'p2'},[1,size(Subregions,2)]),Subregions)==1);

ReMLFRef=zeros(size(FolderPaths.DataFolders,2),size(Subregions,2));ReMLFOLS=zeros(size(FolderPaths.DataFolders,2),size(Subregions,2));ReMLFBinWeight=zeros(size(FolderPaths.DataFolders,2),size(Subregions,2),size(ManualExclusionIndx,2));
maxReMLF=zeros(size(FolderPaths.DataFolders,2),size(Subregions,2));maxReMLFIndx=zeros(size(FolderPaths.DataFolders,2),size(Subregions,2));
ReMLF=zeros(size(FolderPaths.DataFolders,2),size(Subregions,2),size(lambda,2));

AnalDir=reshape(FolderPaths.AnalFolders,[size(lambda,2) size(Subregions,2)]);

for datactr=1:size(FolderPaths.DataFolders,2)
    for regctr=1:size(Subregions,2)
        for lambdactr=1:size(lambda,2)
            eval(['load ' fullfile(RootPath,FolderPaths.DataFolders{datactr},AnalDir{lambdactr,regctr},'SPM.mat')])
            ReMLF(datactr,regctr,lambdactr)=SPM.AL.ReML.F;
        end
        if exist('OLSIndx','var')
            ReMLF(datactr,regctr,:)= ReMLF(datactr,regctr,:) - ReMLF(datactr,regctr,OLSIndx);
        end
        [maxReMLF(datactr,regctr),maxReMLFIndx(datactr,regctr)]=max(ReMLF(datactr,regctr,:));
        ReMLFRef(datactr,regctr)=ReMLF(datactr,regctr,LambdaRefIndx);
        ReMLFOLS(datactr,regctr)=ReMLF(datactr,regctr,OLSIndx);
        ReMLFBinWeight(datactr,regctr,:)=ReMLF(datactr,regctr,ManualExclusionIndx);
    end
    
    if ~isempty(GMIndx) && ~isempty(WMIndx)
        figure
        bar([squeeze(ReMLF(datactr,GMIndx,NotOLS)),squeeze(ReMLF(datactr,WMIndx,NotOLS))]')
        myLegend=cell(size(NotOLS,2),1);
        WLSctr=0;
        for ctrtest=1:size(NotOLS,2)
            if ismember(NotOLS(ctrtest),ManualExclusionIndx)
                ExclFract = regexp(lambda{NotOLS(ctrtest)},'\d*','Match');
                ExclFract = str2num(ExclFract{1});
                myLegend{ctrtest}=['Exclusion - ' num2str(ExclFract) '%'];
            else
                WLSctr=WLSctr+1;
                if size(lambda{NotOLS(ctrtest)},2)==1
                    myLegend{ctrtest}=['WLS \alpha = ' num2str(lambda{NotOLS(ctrtest)}(end))];
                else
                    myLegend{ctrtest}=['WLS \alpha = ' num2str(min(lambda{NotOLS(ctrtest)})) ' to ' num2str(max(lambda{NotOLS(ctrtest)}))];
                end
            end
        end
        legend(myLegend)
        SavePath=fullfile(RootPath,'ReMLAnal',FolderPaths.DataFolders{datactr});
        if ~exist(SavePath,'dir')
            mkdir(SavePath)
        end
        saveas(gcf,fullfile(SavePath,'ReMLgain_histo'),'fig')
    end
end

Indx=[];
for ctr=1:size(Subregions,2)
    if ismember(Subregions{ctr},Params.BrainRegions.WMregions)|strcmp(Subregions{ctr},'p1')|strcmp(Subregions{ctr},'p2')
        Subregions{ctr}=0;
        Indx=cat(2,Indx,ctr);
    end
end
Subregions(Indx)=[];maxReMLF(:,Indx)=[];maxReMLFIndx(:,Indx)=[];ReMLFRef(:,Indx)=[];ReMLFOLS(:,Indx)=[];

if ~isempty(Subregions)
    NMatlas=spm_read_vols(spm_vol(spm_select('FPList',Params.NMDir,'^label.*.nii$')));
    for datactr=1:size(FolderPaths.DataFolders,2)
        ReMLFMaxmap=zeros(size(NMatlas));ReMLFMaxIndxmap=zeros(size(NMatlas));
        ReMLFRefmap=zeros(size(NMatlas));ReMLFOLSmap=zeros(size(NMatlas));
        for regctr=1:size(Subregions,2)
            Indx=find(NMatlas==Subregions{regctr}(1)|NMatlas==Subregions{regctr}(2));
            ReMLFMaxmap(Indx)=maxReMLF(datactr,regctr);ReMLFMaxIndxmap(Indx)=maxReMLFIndx(datactr,regctr);
            ReMLFRefmap(Indx)=ReMLFRef(datactr,regctr);
            ReMLFOLSmap(Indx)=ReMLFOLS(datactr,regctr);
        end
        Vsave=spm_vol(spm_select('FPList',fullfile(RootPath,FolderPaths.DataFolders{datactr},AnalDir{1,1}),'^spmF_0001.nii'));
        
        SavePath=fullfile(RootPath,'ReMLAnal',FolderPaths.DataFolders{datactr});
        if ~exist(SavePath,'dir')
            mkdir(SavePath)
        end
        Vsave.fname=fullfile(SavePath,spm_str_manip(Vsave.fname,'t'));
        MaxToRefmap=ReMLFRefmap-ReMLFMaxmap;MaxToRefmap(find(ReMLFMaxmap==0))=-100;
        Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'ReMLFMax.nii');
        spm_write_vol(Vsave,ReMLFMaxmap);%Map of the maximum Free Energy values, across all noise models
        Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'ReMLFMaxIndx.nii');
        spm_write_vol(Vsave,ReMLFMaxIndxmap);%Map of the model index yielding maximum Free Energy values
        Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'ReMLFRef.nii');
        spm_write_vol(Vsave,ReMLFRefmap);%Map of the Free Energy values for the reference noise models
        Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'ReMLFOLS.nii');
        spm_write_vol(Vsave,ReMLFOLSmap);%Map of the maximum Free Energy values, for OLS analyses
        Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'ReMLFMaxtoRef.nii');
        spm_write_vol(Vsave,MaxToRefmap);%Map of the Free Energy change from the maximum model to the reference
        Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'ReMLFOLStoRef.nii');
        spm_write_vol(Vsave,ReMLFRefmap-ReMLFOLSmap);%Map of the Free Energy change from the OLS analyses to the reference noise model
    end
end
end