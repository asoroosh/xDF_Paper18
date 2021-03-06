
clear;

%%%%%%%%%%%%%%%%%%%%%%%%
% This script put the ACI values on the Atlases used in the paper. 
% You need the ACI previousely calculated: e.g.
% 'R/HCP100UR/' a{1} '_' g{1} '_FPP_CorrLeng_ROIWise.mat'
%
%%%%%%%%%%%%%%%%%%%%%%%%

Atlas={'ICA200','Power','Yeo','Gordon'};
GSR={'GSR'};

% msk_rt='/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/Atlas/';
% msk_rt has to point into the directory that contains all the Atlases used
% in the paper.

for a=Atlas
    for g=GSR
    
        load(['R/HCP100UR/' a{1} '_' g{1} '_FPP_CorrLeng_ROIWise.mat']);
        
        if strcmp(g{1},'NoGSR')
            cl=mean(ACL,2); %aal2 gsr
        elseif strcmp(g{1},'GSR')
            cl=mean(ACL,2); %aal2 gsr
        end

        if strcmp(a{1},'Power')
            niiAt=load_untouch_nii([msk_rt 'Power_ROIs/Power2011_MNI2mm_AgrTrmd.nii.gz']);
        elseif strcmp(a{1},'Yeo')
            niiAt=load_untouch_nii([msk_rt 'Yeo2011_17Networks_FSL_MNI152_2mm.nii.gz']);
        elseif strcmp(a{1},'Gordon')
            niiAt=load_untouch_nii([msk_rt 'Parcels_MNI_222.nii']);  
        %elseif strcmp(a{1},'ICA200')
        %    niiAt=load_untouch_nii(['/Users/sorooshafyouni/Desktop/HCP/PTN/HCP_PTN1200/groupICA/groupICA_3T_HCP1200_MSMAll_d200.ica/ICA200_thr_z_25.nii.gz']);   
        end

        mimg=niiAt.img;
        for r=1:numel(cl)
            [xx,yy,zz]=ind2sub(size(mimg),find(mimg==r));
            for vx=1:length(xx)      
                mimg(xx(vx),yy(vx),zz(vx))=cl(r);
            end
        end

        niiAt.img=mimg;
        save_untouch_nii(niiAt,['R/CorrLength_on_' a{1} '_' g{1} '.nii'])

        clear niiAt mimg cl xx yy zz
    end
end

