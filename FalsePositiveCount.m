function FalsePositiveCount(FolderPaths,DataStr,Thr1,Thr2)
% Tracks the occurence of significant results across multiple repetitions of a statistical analysis.
% INPUTS:
%   - FolderPaths: structure containing the paths to the analysis repetitions folder - computed in PrepAnalysis;
% 	- DataStr: path to the analysis folder, assumed to be within the repetition folders
%   - Thr1/Thr2: significance thresholds for voxel and cluster-level analyses
% OUTPUTS: saved to disk. In 'FalsePosAnalysis' folder. 
%     - FalsePosRuns_clus.mat. Number of repetitions with false positives at the cluster level. 
%     - FalsePosRuns_vxl.mat. Number of repetitions with false positives at the voxel level.
%
% Uses the cp_cluster_Pthresh.m routine written by Christophe Phillips, Cyclotron Research Centre, University of Liege, Belgium 
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

for ctr=1:size(FolderPaths.DataFolders,2)
    for datactr=1:size(FolderPaths.AnalFolders,2)
        SavePath=fullfile(spm_str_manip(FolderPaths.CohortPaths{1},'h'),'FalsePosAnalysis',FolderPaths.DataFolders{ctr},FolderPaths.AnalFolders{datactr});
        if ~exist(SavePath,'dir')
            mkdir(SavePath)
        end
        FalsePosRuns_vxl=0;Nbmissing=0;NbClusters=0;FalsePosRuns_clus=0;
        for repeatctr=1:size(FolderPaths.CohortPaths,2)
            CurrentPath=fullfile(FolderPaths.CohortPaths{repeatctr},FolderPaths.DataFolders{ctr},FolderPaths.AnalFolders{datactr});
            if exist(CurrentPath,'dir')==7 && ~isempty(spm_select('FPList',CurrentPath,DataStr))
                StatFiles=spm_select('FPList',CurrentPath,DataStr);Statmap=abs(spm_read_vols(spm_vol(StatFiles(1,:))));
                if ~isempty(find(Statmap>Thr1(datactr)))%voxel-level
                    FalsePosRuns_vxl=FalsePosRuns_vxl+1;
                    if FalsePosRuns_vxl==1
                        FalsePosMap=zeros(size(Statmap));
                        Vsave=spm_vol(spm_select('FPList',CurrentPath,DataStr));Vsave=Vsave(1); 
                    end
                    FalsePosMap(find(Statmap>Thr1(datactr)))=FalsePosMap(find(Statmap>Thr1(datactr)))+1;
                end
                if ~isempty(find(Statmap>Thr2(datactr)))%cluster-level
                    testctr=1;
                    while testctr~=size(StatFiles,1)+1
                        xSPM.Ic=testctr;xSPM.u=Thr2(datactr);xSPM.Im=[];
                        xSPM.k = 1;xSPM.thresDesc='none';
                        xSPM.swd= CurrentPath;
                        [SPM,xSPM] = spm_getSPM(xSPM);% A = spm_clusters(xSPM.XYZ);%unique(A)
                        cp_clust_size=cp_cluster_Pthresh(xSPM,0.05);
                        clear xSPM;
                        xSPM.Ic=testctr;xSPM.u=Thr2(datactr);xSPM.Im=[];
                        xSPM.k = cp_clust_size;xSPM.thresDesc='none';
                        xSPM.swd= CurrentPath;
                        [SPM,xSPM] = spm_getSPM(xSPM);
                        A = spm_clusters(xSPM.XYZ);%unique(A)
                        if ~isempty(A)
                            testctr=size(StatFiles,1)+1;
                            NbClusters=NbClusters+size(unique(A),2);%max(A);sum(A==3);max(xSPM.Z(A==1));
                            FalsePosRuns_clus=FalsePosRuns_clus+1;
                        else
                            testctr=testctr+1;
                        end
                    end
                end
            else
                Nbmissing=Nbmissing+1;
            end
        end
        if exist('FalsePosMap','var')
            Vsave.fname=fullfile(SavePath,'FalsePositives.nii');
            spm_write_vol(Vsave,FalsePosMap);
            FalsePosSize=sum(FalsePosMap(:));
        end
        save(fullfile(SavePath,'NbClusters'),'NbClusters', '-v7.3')
        save(fullfile(SavePath,'FalsePosRuns_vxl'),'FalsePosRuns_vxl', '-v7.3')
        save(fullfile(SavePath,'FalsePosRuns_clus'),'FalsePosRuns_clus', '-v7.3')
        if exist('FalsePosSize','var')
            save(fullfile(SavePath,'FalsePosSize'),'FalsePosSize', '-v7.3')
        end
    end
end

end