clear

addpath /Users/sorooshafyouni/Home/Atlas

PCC = load_untouch_nii(['/Users/sorooshafyouni/Home/Eff_Con/RSN10/Comps/PPC_mask_RSN10.nii.gz']);
aimg = PCC.img;

simg = load_untouch_nii(['/Users/sorooshafyouni/Desktop/HCP/135932/135932_RS/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz']);
simg = simg.img;

[mTS,ROI,SanCheck]=tsExtraction(simg,aimg);

[Y,Stat]=CleanNIFTI(simg,'demean');

ARC = 3;
[rSBfc,pSBfc] = corr(mTS(ARC:end),Y(:,1:end-(ARC-1))');
%fdr_bh(pSBfc)

I0 = Stat.OrigDim(1);
T = Stat.OrigDim(2);

[Xd,Yd,Zd,Td] = size(simg);
idx        = 1:I0;
idx(Stat.Removables) = [];

Var1_tmp = zeros(I0,1);
Var1_tmp(idx) = rSBfc;
rSBfc = reshape(Var1_tmp,[Xd Yd Zd]);

Var1_tmp = zeros(I0,1);
Var1_tmp(idx) = pSBfc;
pSBfc = reshape(Var1_tmp,[Xd Yd Zd]);

WindSize = 3;
pSBfc = imgaussfilt3(pSBfc,1,'FilterSize',[WindSize WindSize WindSize]);

SBfc = rSBfc.*fdr_bh(pSBfc);

nii_tmp = make_nii(SBfc,[2,2,2],[0,0,0],64,['3D image of Seed-based FC']);
save_nii(nii_tmp,['SBfc_AR' num2str(ARC) '.nii.gz'])
