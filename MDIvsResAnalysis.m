function MDIvsResAnalysis(QUIQI,FolderPaths)
% Analysis of image noise vs motion degradation index.
% INPUTS:
%     - QUIQI: structure containing all information used for analysis. Computed in PrepAnalysis.m.
%     - FolderPaths: path to the demographic information structure. Computed in PrepAnalysis.m.
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

eval(['load ' fullfile(FolderPaths.CohortPaths{:},'Subject_Details.mat')]);
MDIVals=[];
for ctr=1:size(Subject_Details,2)
    MDIVals=cat(1,MDIVals,[Subject_Details(ctr).QA.SDR2s.MTw Subject_Details(ctr).QA.SDR2s.PDw Subject_Details(ctr).QA.SDR2s.T1w]);
end

for datactr=1:size(QUIQI,2)    
    MDIVals4Anal{datactr}=MDIVals(:,QUIQI(datactr).SDR2sIndx);
    
    CurrentPath=fullfile(QUIQI(datactr).CohortPath,QUIQI(datactr).AnalDir);
    P=spm_select('FPList',CurrentPath,'^ExplicitMask_.*.nii$');
    Vsave=spm_vol(P);
    Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'ReMLMask.nii');
    load(fullfile(CurrentPath,'SPM.mat'));%reads-in ReMLMask from SPM matrix and saves it as nifti for residual analysis below
    spm_write_vol(Vsave,SPM.AL.ReML.Mask);
end

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end
parpool;
parfor datactr=1:size(QUIQI,2)
% for datactr=1:size(ROI,1)
    CurrentPath=fullfile(QUIQI(datactr).CohortPath,QUIQI(datactr).AnalDir);
    ResFiles=spm_select('FPList',CurrentPath,'^Res_.*.nii$');
    Vsave=spm_vol(ResFiles(1,:));
    
    ExplicitMask=spm_read_vols(spm_vol(spm_select('FPList',CurrentPath,'^ExplicitMask_.*.nii$')));% for KS analysis
    ReMLMask=spm_read_vols(spm_vol(spm_select('FPList',CurrentPath,'ReMLMask.nii'))); %for ResvsMDI analysis
    ExplicitMaskIndx=find(ExplicitMask~=0);ReMLMaskIndx=find(ReMLMask~=0);
    
    ResidVar=zeros(size(MDIVals,1),1);Residuals=zeros(size(ExplicitMaskIndx,1),size(MDIVals,1));
    for Subjctr=1:size(MDIVals,1)% reads-in individual residual maps and estimates variance
        Subjctr
        tempRes=spm_read_vols(spm_vol(ResFiles(Subjctr,:)));
        ResidVar(Subjctr)=var(tempRes(ReMLMaskIndx),'omitnan');% For residual analysis. Considers the voxels used by ReML for noise modelling 
        Residuals(:,Subjctr)=tempRes(ExplicitMaskIndx);% For KS analysis below. Considers all voxels of the explicit mask used for analysis
    end
    SavePath=fullfile(spm_str_manip(Vsave.fname,'h'),'ResidualAnalysis');
    if ~exist(SavePath,'dir')
        mkdir(SavePath)
    end
    Vsave.fname=fullfile(SavePath,spm_str_manip(Vsave.fname,'t'));
    parsave(fullfile(SavePath,'ResidVar'),ResidVar)
    
    if ~isempty(strfind(QUIQI(datactr).AnalDir,'R2s'))
        MDIvsResFit(SavePath,MDIVals4Anal{datactr},ResidVar*1e6)
        parsave(fullfile(SavePath,'MDIVals'),MDIVals4Anal{datactr})
    elseif ~isempty(strfind(QUIQI(datactr).AnalDir,'R1'))
        % %     bi-variate
        %     DM=cat(2,ones(size(MDIVals,1),1),MDIVals(:,2:3));
        %     beta=mvregress(DM,ResidVar(:,2));
        %     SSresid = sum((ResidVar(:,2) - DM*beta).^2);
        %     SStotal = (length(ResidVar(:,2))-1) * var(ResidVar(:,2));
        %     Rsq = 1 - SSresid/SStotal;
        %
        %     figure
        %     plot(ResidVar(:,2),'.')
        %     hold on
        %     plot(DM*beta,'r.')
        %     title(['(\alpha_1, \alpha_2 ) = (' num2str(round(beta(2)*1e1)/1e1) ', ' num2str(round(beta(3)*1e1)/1e1) '); R^2 = ' num2str(round(Rsq*1e2)/1e2)])
        %     saveas(gcf, fullfile(LocalPath,char(DataDir),'Residuals_bivariate'), 'fig');
        %     close(gcf)
        %
        %     univariate
        
        MDIvsResFit(NMatlas,ExplicitMask,TempAtlasVals,Vsave,SavePath,sum(MDIVals4Anal{datactr},2),ResidVar,LocalResidVar)
        
        %             [P,Rsq]=myPolyFit(sqrt(sum(MDIVals(:,2:3).^2,2)),ResidVar,2);
        %             plotLinFit(sqrt(sum(MDIVals(:,2:3).^2,2)),ResidVar',P,Rsq,'\surd (SoS(SDR2s)) (s^-^1)','Residuals (Var)')
        %             saveas(gcf, fullfile(SavePath,'Residuals_SDR2sSoS'), 'fig');close(gcf)
    elseif ~isempty(strfind(QUIQI(datactr).AnalDir,'MT'))
        %     univariate
        [P,Rsq]=myPolyFit(sum(tempSDR2sVals,2),ResidVar(:,1),2);
        plotLinFit(sum(MDIVals,2),ResidVar(:,1)',P,Rsq,'\Sigma SDR2s (s^-^1)','Residuals (Var)')
        saveas(gcf, fullfile(SavePath,'Residuals_SDR2sSum'), 'fig');close(gcf)
        [P,Rsq]=myPolyFit(sqrt(sum(tempSDR2sVals.^2,2)),ResidVar(:,1),2);
        plotLinFit(sqrt(sum(tempSDR2sVals.^2,2)),ResidVar(:,1)',P,Rsq,'\surd (SoS(SDR2s)) (s^-^1)','Residuals (Var)')
        saveas(gcf, fullfile(SavePath,'Residuals_SDR2sSoS'), 'fig');close(gcf)
    end

%     KS tests
    HistResH=zeros(size(Residuals,1),1);HistResP=zeros(size(Residuals,1),1);
    for ctr=1:size(Residuals,1)%1e5%
        Signal=Residuals(ctr,:)';
        if isempty(find(isnan(Signal)))&&isempty(find(isinf(Signal)))
            [HistResH(ctr),HistResP(ctr)]=kstest((Signal-mean(Signal))/std(Signal,[],1));
        end
    end
    Hmap=zeros(size(spm_read_vols(spm_vol(ResFiles(1,:)))));
    Hmap(ExplicitMaskIndx)=HistResH;    
    Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'HistResH.nii');spm_write_vol(Vsave,Hmap);

    Pmap=zeros(size(spm_read_vols(spm_vol(ResFiles(1,:)))));
    Pmap(ExplicitMaskIndx)=HistResP;    
    Vsave.fname=fullfile(spm_str_manip(Vsave.fname,'h'),'HistResP.nii');spm_write_vol(Vsave,Pmap); 
    
    for ctr=1:size(ResFiles,1)%delete residual individual residual maps to save disk space
        delete(deblank(ResFiles(ctr,:)));
    end
    
end

end

function MDIvsResFit(SavePath,MDI,Res)
Params=GetParams();

if ~exist(SavePath,'dir')
    mkdir(SavePath)
end

[P,Rsq,yfit]=myPolyFit(MDI,Res,Params.MDIvsResOrder,'NonNeg');
% [P,Rsq,yfit]=myPolyFit(MDI,Res,Params.MDIvsResOrder,'Free');
plotLinFit(MDI,Res',yfit,P,Params.MDIvsResOrder,Rsq,'MDI (s^-^1)','Residuals (Var)')
saveas(gcf, fullfile(SavePath,'Residuals'), 'fig');
close(gcf)
FittingEstimates.P=P;FittingEstimates.Rsq=Rsq;
save(fullfile(SavePath,'FitEstimates.mat'),'FittingEstimates', '-v7.3')

end
function parsave(Path,Var)
save(Path,'Var')
end


