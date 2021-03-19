function [Subject_Details_subset]=Select_cohort_subset(Subject_Details,AnalParams)
% From a full cohort of data, selects a subset of data with the same number of datasets per age bin
% INPUTS:
%     - Subject_Details: structure of demographic information for the analysis cohort 
%     - AnalParams: structure containg the details of the binning of the
%     input cohort according to age. Set in GetAnalParams.m
% OUTPUTS:
%     -  Subject_Details_subset: subset of data with AnalParams.AgeBinning.Nsamples datasets per age bin
%
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

[BinSize,BinIndx]=histc([Subject_Details.Age],linspace(AnalParams.AgeBinning.Agemin,AnalParams.AgeBinning.Agemax,AnalParams.AgeBinning.BinNb));

Subject_Details_subset=[];
for ctr=1:AnalParams.AgeBinning.BinNb
    Indx=find(BinIndx==ctr);
    temp=randperm(size(Indx,2));temp=temp(1:min(AnalParams.AgeBinning.Nsamples,BinSize(ctr)));
    Subject_Details_subset=cat(2,Subject_Details_subset,...
        Subject_Details(Indx(temp)));
end
end

