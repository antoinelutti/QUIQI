function Subregions=getROIpairs(TissueType)
% Computes pairs of regional atlas labels extracted from the left & right hemispheres.
% INPUTS:
%     - TissueType: tissue class of interest for the extraction of the atlas labels. Possible values: 'WM' (white
% matter, 'GM' (grey matter) or 'GMnoB0' (grey matter after exclusion of the regions affected by B0-inhomogeneities, defined in GetParams)

% OUTPUTS:
%     - Subregions: cell array of pairs of regional atlas labels (left/right hemisphere). 
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

Params=GetParams;
NMatlas=spm_read_vols(spm_vol(spm_select('FPList',Params.NMDir,'^label.*.nii$')));
NMvals=unique(NMatlas);
NMvals=NMvals(find(~ismember(NMvals,Params.BrainRegions.NotGMorWM) & ~ismember(NMvals,[23 30 35 69 71 72 73])));
if nargin>0%selects one tissue type only
    if strcmp(TissueType,'GM')
        NMvals=NMvals(find(ismember(NMvals,Params.BrainRegions.GMregions)));
    elseif strcmp(TissueType,'GMnoB0')
        NMvals = NMvals(ismember(NMvals,Params.BrainRegions.GMregions)&~ismember(NMvals,Params.BrainRegions.B0regions));
    elseif strcmp(TissueType,'WM')
        NMvals=NMvals(find(ismember(NMvals,Params.BrainRegions.WMregions)));
    end
end
Evenvals=NMvals(1:2:end);Oddvals=NMvals(2:2:end);
Subregions=cell(1,size(Evenvals,1));
for ctr=1:size(Evenvals,1)
    Subregions(ctr)={cat(1,Evenvals(ctr),Oddvals(ctr))};
end

end