function DM=MakeDM(Subject_Details,MDI,AnalType)
% Computes a design  matrix DM used in subsequent analyses (RunAnalysis.m). 
% INPUTS:
%     - Subject_Details: structure of demographic information for the analysis cohort e.g. age,
% gender, brain volume,... 
%     - MDI: Motion Degradation Index values for each dataset of the analysis cohort. Inserted in the design matrix for
% 'MotionBias' analyses, for each power value of Params.MotionRegPowers 
%     - AnalType: string of the analysis type 

% OUTPUTS:
%     - DM.mat: design matrix data
%     - DM.size: size. Used for the definition of the F-contrast in the
%     ageing analyses 
%     - DM.text: text description
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland


Params=GetParams;
Covariates=[Subject_Details.confound];
if strcmp(AnalType,'MotionBias')
    MDI = MDI-mean(MDI,1);%mean centering and power components is done here to allow for analysis on a subset of the cohort  
    DM.mat=MDI.^Params.MotionRegPowers(1);
    for ctr=2:size(Params.MotionRegPowers,2)
        DM.mat=cat(2,DM.mat,MDI.^Params.MotionRegPowers(ctr));    
    end
    DM.size=size(DM.mat,2);
    for ctr=1:size(DM.mat,2)
        DM.text{ctr}=num2str(ctr);
    end
elseif strcmp(AnalType,'GroupComparison')
    DM.mat=[];DM.size=0;
    DM.mat=cat(2,DM.mat,[Subject_Details(:).Age]'-mean([Subject_Details(:).Age],2),...
        ([Subject_Details(:).Age]'-mean([Subject_Details(:).Age],2)).^2);
    DM.text={'Age','Agesq'};
else
    DM.mat=cat(2,[Subject_Details(:).Age]',([Subject_Details(:).Age]-mean([Subject_Details.Age])).^2');
    DM.text={'Age','Agesq'};
    DM.size=size(DM.text,2);
end
DM.mat=cat(2,DM.mat,[Covariates.Gender]');DM.text(size(DM.mat,2))={'Gender'};

DM.desc=AnalType;
end
