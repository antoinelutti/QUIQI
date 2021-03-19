function Name=RegionLabeltoName(Label)
% Creates the name of a folder that will contain the analysis results  
% INPUTS:
%     -  Label: region of interest in the analysis. Initialized in RunQUIQI.m.
% OUTPUTS:
%     -  Name: cell array of strings.
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland
Name=cell(size(Label,2),1);
for ctr=1:size(Label,2)
    if ischar(Label{ctr})
        Name{ctr}=Label{ctr};
    else
        Delimiter=cell(1,size(Label{ctr},1)-1);
        for ctr2=1:size(Delimiter,2)
            Delimiter{ctr2}='_';
        end
        Name{ctr,1}=strjoin(cellstr(num2str(Label{ctr})),Delimiter);
    end
end

end