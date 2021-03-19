function Params=GetParams()
% Returns a structure containing details of the analysis to be conducted. 
% 
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

Params.HomeDir='F:\alutti\QUIQIdata';% The folders containing the analysis results will be created at this location
Params.DataDir=Params.HomeDir;% Location of the data folders
Params.NMDir=fullfile(Params.HomeDir,'Neuromorphometrics');% Location of the Neuromorphometrics folder
Params.AnalysisFolder='Analysis';%Naming of the analysis folder
Params.DataSubDir='';%Data location within the main folder of each dataset.

Params.MDIvsResOrder=4;%order of the polynomial fit of var(res)vsMDI

Params.ReMLAnal.RefLambda=3;%power of the reference noise model for analysis of the Free Energy change
Params.MotionRegPowers=[1 2 3 4];% powers of the MDI inserted in design matrix for MotionBias analyses

% labelling of grey and white matter regions
NMatlas=spm_read_vols(spm_vol(spm_select('FPList',Params.NMDir,'^label.*.nii$')));
AtlasVals=unique(NMatlas);
Params.BrainRegions.NotGMorWM=[0 4 11 46 49 50 51 52 61 62];
Params.BrainRegions.WMregions=[35 40 41 44 45];
Params.BrainRegions.GMregions=AtlasVals;Params.BrainRegions.GMregions(ismember(Params.BrainRegions.GMregions,Params.BrainRegions.WMregions)...
    |ismember(Params.BrainRegions.GMregions,Params.BrainRegions.NotGMorWM))=[];
Params.BrainRegions.B0regions=[102 103 104 105 120 121 122 123 124 125 132 133 140 141 146 147 154 155 178 179 186 187 202 203 ...%OFC
31 32 47 48 116 117 170 171];%amygdala/temporal lobe/hippocampus

end