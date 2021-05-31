function AnalParams=GetAnalParams(AnalType)
% Retrieves parameters specific to each analysis type.
% INPUTS: 
%     - AnalType: string label of the analysis to be conducted.
% OUTPUTS:
%     - AnalParams: structure of analysis parameters. 
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

% Default analysis settings are set for a standard 'Full' analysis of the whole cohort.
% For other analysis types, the relevant parameters are set specifically.
AnalParams.SaveResiduals=0;AnalParams.Subcohort=0;AnalParams.Exclusion=0;
AnalParams.AgeBinning.Nsamples=Inf;AnalParams.NRepeats=1;
AnalParams.AgeBinning.Agemin=10;AnalParams.AgeBinning.Agemax=100;AnalParams.AgeBinning.BinNb=19;
AnalParams.ShuffleAge=0;AnalParams.ShuffleData=0;AnalParams.NGroup1=0;

if strcmp(AnalType,'Exclusion')
    AnalParams.Exclusion=[0 3 7 13 20 30];
    AnalParams.NRepeats=size(AnalParams.Exclusion,2);
elseif strcmp(AnalType,'Residuals')
    AnalParams.SaveResiduals=1;
elseif strcmp(AnalType,'Reproducibility')
    AnalParams.Subcohort=1;AnalParams.AgeBinning.Nsamples=200;AnalParams.NRepeats=100;
elseif strcmp(AnalType,'Specificity')
    AnalParams.AgeBinning.Nsamples=10;AnalParams.NRepeats=1e1;
    AnalParams.Subcohort=1;AnalParams.ShuffleAge=1;
elseif strcmp(AnalType,'MotionBias')
    AnalParams.SaveResiduals=1;
    AnalParams.Subcohort=1;AnalParams.AgeBinning.Nsamples=200;
    AnalParams.AgeBinning.Agemin=56;AnalParams.AgeBinning.Agemax=58;AnalParams.AgeBinning.BinNb=2;
elseif strcmp(AnalType,'GroupComparison')
    AnalParams.Subcohort=1;AnalParams.ShuffleData=1;AnalParams.NGroup1=60;
    AnalParams.AgeBinning.Nsamples=200;AnalParams.NRepeats=1e1;
    AnalParams.AgeBinning.Agemin=56;AnalParams.AgeBinning.Agemax=58;AnalParams.AgeBinning.BinNb=2;
end
end