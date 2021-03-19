function [QUIQI,Folders]=PrepAnalysis(Subject_Details_All,lambda,DataTypes,Subregions,AnalType)
% Compiles into QUIQI and Folders structures all information required for
% subsequent analysis. Also, edits cohort details if necessary (e.g. for Exclusion or Specificity analyses)
% and computes an explicit mask for image analysis if it doesn't exist (PrepareMasksAndCohort function)

% INPUTS: 
%     - Subject_Details_All: Subject demographic information structure.
%     - DataTypes: type of input data for the analysis. Initialized in RunQUIQI.m.
%     - Subregions: region of interest of the analysis. Initialized in RunQUIQI.m.
%     - AnalType: analysis type. Initialized in RunQUIQI.m.

% OUTPUTS: 
% QUIQI.TissueType - tissue type of interest in the analysis
% QUIQI.SDR2sIndx - index of the MDI of interest in the demographic information structure, which
% contains three entries. Matches the type of data to be analysed i.e. 1/2/3 for
% MT/PD/T1-weighted data
% QUIQI.MotionReg - MDI values extracted from the demographic information structure. From the three entries, the extracted values are taken from the value of QUIQI.SDR2sIndx. Used for the MotionBias analyses.
% QUIQI.ReML - set of basis functions for ReML
% QUIQI.ROI - region of interest
% QUIQI.AnalDir - analysis output folder
% QUIQI.InputData - Input data string 

%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

Params=GetParams;
AnalParams=GetAnalParams(AnalType);
RegionStr=RegionLabeltoName(Subregions);
NbRepeats=AnalParams.NRepeats;

RootPath=fullfile(Params.HomeDir,[Params.AnalysisFolder '_' AnalType]);
if ~exist(RootPath,'dir')
    mkdir(RootPath)
end
save(fullfile(RootPath,'AnalParams'),'AnalParams', '-v7.3')
%     COHORT DEFINITION
CreateNewMask=zeros(NbRepeats,1);CohortPath=cell(1,NbRepeats);
for subsetctr=1:NbRepeats
    CohortPath{subsetctr}=fullfile(RootPath,num2str(subsetctr));
    if ~exist(CohortPath{subsetctr},'dir')
        mkdir(CohortPath{subsetctr})
    end
    if isempty(spm_select('FPList',CohortPath{subsetctr},'^Subject_Details.*.mat$'))
        if ~isinf(AnalParams.AgeBinning.Nsamples)
            Subject_Details=Select_cohort_subset(Subject_Details_All,AnalParams);
        else
            Subject_Details=Subject_Details_All;
        end
        if strcmp(AnalType,'Exclusion') 
            if AnalParams.Exclusion(subsetctr)~=0%Applicable to PDw only (for simplicity)
                PDSDR2sVals=zeros(size(Subject_Details,2),1);
                for subjctr=1:size(Subject_Details,2)
                    PDSDR2sVals(subjctr)=Subject_Details(subjctr).QA.SDR2s.PDw;
                end
                [~,B]=sort(PDSDR2sVals);
                Indx=B(end-round(size(PDSDR2sVals,1)*AnalParams.Exclusion(subsetctr)/100)+1:end);
                Subject_Details(Indx)=[];
            end
        end
        if AnalParams.ShuffleAge
            Age=[Subject_Details(:).Age];
            ShuffledAge=Age(randperm(size(Subject_Details,2)));
            for ctr=1:size(Subject_Details,2)
                Subject_Details(ctr).Age=ShuffledAge(ctr);
            end
        elseif AnalParams.ShuffleData
            Subject_Details=Subject_Details(randperm(size(Subject_Details,2)));
        end
        CreateNewMask(subsetctr)=1;
        save(fullfile(CohortPath{subsetctr},'Subject_Details'),'Subject_Details', '-v7.3')
    end
end
%     END OF COHORT DEFINITION

%     MASK CREATION
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool(12);
parfor subsetctr=1:NbRepeats
    PrepareMasksAndCohort(CohortPath{subsetctr},CreateNewMask(subsetctr),Subregions)
end

Folders.CohortPaths=CohortPath;
Folders.AnalFolders={};Folders.DataFolders={};
%     BUILD QUIQI structures
ctr1d=0;
for subsetctr=1:NbRepeats
    eval(['load ' spm_select('FPList',CohortPath{subsetctr},'^Subject_Details.*.mat$')]);
    SDR2sVals=zeros(size(Subject_Details,2),3);
    for subjctr=1:size(Subject_Details,2)
        SDR2sVals(subjctr,:)=[Subject_Details(subjctr).QA.SDR2s.MTw Subject_Details(subjctr).QA.SDR2s.PDw Subject_Details(subjctr).QA.SDR2s.T1w];
    end
    for ctr1=1:size(DataTypes,2)
        for ctr=1:size(Subregions,2)
            for ctr2=1:size(lambda,2)
                ctr1d=ctr1d+1;
                QUIQI(ctr1d).CohortPath=CohortPath{subsetctr};
                if ischar(Subregions{ctr})
                    QUIQI(ctr1d).TissueType=Subregions(ctr);
                elseif isnumeric(Subregions{ctr})
                    if ismember(Subregions{ctr},Params.BrainRegions.GMregions)==1
                        QUIQI(ctr1d).TissueType={'p1'};
                    elseif ismember(Subregions{ctr},Params.BrainRegions.WMregions)==1
                        QUIQI(ctr1d).TissueType={'p2'};
                    end
                end
                if strcmp(DataTypes{ctr1},'MTw_R2s')
                    QUIQI(ctr1d).SDR2sIndx=1;
                elseif strcmp(DataTypes{ctr1},'PDw_R2s')
                    QUIQI(ctr1d).SDR2sIndx=2;
                elseif strcmp(DataTypes{ctr1},'T1w_R2s')
                    QUIQI(ctr1d).SDR2sIndx=3;
                elseif strcmp(DataTypes{ctr1},'PDT1_R2s')
                    QUIQI(ctr1d).SDR2sIndx=[2 3];
                elseif strcmp(DataTypes{ctr1},'R2s_OLS')
                    QUIQI(ctr1d).SDR2sIndx=[1 2 3];
                elseif strcmp(DataTypes{ctr1},'R1')
                    QUIQI(ctr1d).SDR2sIndx=[2 3];
                elseif strcmp(DataTypes{ctr1},'MT')
                    QUIQI(ctr1d).SDR2sIndx=[1 2 3];
                end
                QUIQI(ctr1d).MotionReg={SDR2sVals(:,QUIQI(ctr1d).SDR2sIndx)};
                if isnumeric(lambda{ctr2})
                    TempReML=[];
                    for ctr3=1:size(lambda{ctr2},2)
                        if lambda{ctr2}(ctr3)==0% to avoid duplicates of the identity matrix
                            TempReML=cat(2,TempReML,ones(size(SDR2sVals,1),1));
                        else
                            TempReML=cat(2,TempReML,SDR2sVals(:,QUIQI(ctr1d).SDR2sIndx).^lambda{ctr2}(ctr3));
                        end
                    end
                    QUIQI(ctr1d).ReML=TempReML;
                    LambdaStr=num2str(lambda{ctr2});LambdaStr=LambdaStr(find(isspace(LambdaStr)==0))';FolderName={};
                    for strctr=1:size(LambdaStr,1)
                        FolderName{strctr}=LambdaStr(strctr);
                    end
                    FolderName=['_lambda_' strjoin(FolderName,'_')];
                end
                QUIQI(ctr1d).ROI=RegionStr{ctr};
                QUIQI(ctr1d).AnalDir=fullfile(DataTypes{ctr1}, [char(RegionStr{ctr}) FolderName]);
                QUIQI(ctr1d).InputData=['fin_dart_' char(QUIQI(ctr1d).TissueType) '_ws.*' DataTypes{ctr1} '.nii$'];
                if ctr1==1&&subsetctr==1
                    Folders.AnalFolders(size(Folders.AnalFolders,2)+1)={[char(RegionStr{ctr}) FolderName]};
                end
            end
        end
    end
end
Folders.DataFolders=DataTypes;

end


