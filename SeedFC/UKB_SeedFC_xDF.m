clear

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/DVARS/Aux')) %for CleanNIFTI and the GSRme
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/xDF')) %xDF
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/Utils/Atlases')) %xDF

%%% set up FSL
fsldir = getenv('FSLDIR');
if isempty(fsldir) 
    error('CleanNIFTI_fsl:: I can not find FSL directory!'); 
end
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

%%% Read the NIFTI file

Path2Sub = ['/Users/sorooshafyouni/Desktop/HCP/955465/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz'];

%%% Extract the Seed time series ==========
seed_img = read_avw('seeds/Yeo_LR_PCC.nii.gz'); 
[sub_img,sub_dim,sub_voxelsizes]  = read_avw(Path2Sub);

seed_ts = tsExtraction(sub_img,seed_img);
seed_ts = mean(seed_ts,2);

%seed_ts = seed_ts(:,2);
%%% Reshape and global signal regression (?)

[sub_ts,sub_ts_stat] = CleanNIFTI_fsl(sub_img);
T  = sub_ts_stat.OrigDim(2);
GS = sub_ts_stat.GlobalMeanSignal;

%%% Regress out the global signal from the image and the seed:
disp('== Global signal regression...')
sub_ts  = GSRme(sub_ts,T,'GS',GS');
seed_ts = GSRme(seed_ts,T,'GS',GS');

%%% Do the xDF =============================
disp('==xDF...')
[~,xDF_Stat]=xDF_seedFC(sub_ts',seed_ts',T,'truncate','adaptive','verbose');
Z_xDF = xDF_Stat.z;

%%% Do Pearson's ===========================
disp('==Pearsons correlations...')
R = corr(seed_ts,sub_ts);

if sum(isnan(R)); error('There is something wrong!'); end

Z_Naive = atanh(R);
Z_Naive = Z_Naive.*sqrt(T-3); 


%%% Save the Images ========================

CleanNIFTI_fsl(Z_Naive,'destdir','Ztest.nii.gz','Removables',sub_ts_stat.Removables,'ImgDim',sub_ts_stat.ImgDim,'voxelsize',sub_voxelsizes);
