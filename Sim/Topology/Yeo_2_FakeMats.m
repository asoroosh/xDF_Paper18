clear

j = str2double(getenv('SGE_TASK_ID'));

rng(j)

TVFlag = 'TVOff';

thrng=0.01:0.01:0.50;

addpath(genpath('~/bin/2017_01_15_BCT'))
addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/contest'))

%SubID = 135932;
SubID = 118528;

load(['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FC/Yeo/FPP/HCP_FPP_' num2str(SubID) '_ROIsTS.mat'],'mts')

T  = size(mts,1);
N = size(mts,2);

AC = AC_fft(mts,T);
AC(:,[1 end])=[];

AC = AC(:,1:20); %because you don't expect the full length is real true pure and dignified autocorrelation, right?

% AUTOCORRELATED NETWORKS =============================================================================

    disp(['Net: ' num2str(j)])
    
    for i=1:N
        CovMats(:,:,i)=MakeMeCovMat_sqrtm(AC(i,:),T);
    end
    simts = corrautocorr_sqrtm(zeros(1,N),eye(N),CovMats,T);
    
    clear CovMats
    
    %---MAT ==========================================
    Mat.amat = atanh(corr(simts')).*sqrt(T-3);
    Mat.amat(1:N+1:end) = 0;
    %---

%=====================================================
%=====================================================
% ADJUSTED NETWORKS ==================================
    disp('AC corrections...')
    %---MAT ==========================================
    [~,adjCR_amat] = CRBCF(simts,T);
    Mat.adjCR_amat = adjCR_amat;
    Mat.adjCR_amat(1:N+1:end) = 0;
    
    %---
    %-- ==============================================
    [~,adjAR1_amat] = Bartlett46_fft(simts,T);
    Mat.adjAR1_amat = adjAR1_amat;
    Mat.adjAR1_amat(1:N+1:end) = 0;
    
    %-- ==============================================
    [~,adjCH_amat] = CheltonBCF(simts,T);
    Mat.adjCH_amat = adjCH_amat;
    Mat.adjCH_amat(1:N+1:end) = 0;
    
    % xDF --------------------------------------------------------------------------------------
   disp('xDF corrections...')

    %---MAT ==========================================
    [~,Stat]=xDF(simts,T,TVFlag);
    Mat.adj_amat = Stat.z.rzf; 
    Mat.adj_amat(1:N+1:end) = 0;
    clear Stat
    %---MAT ==========================================
    [~,Stat]=xDF(simts,T,'taper','shrink',TVFlag);
    Mat.adjs1_amat = Stat.z.rzf;
    Mat.adjs1_amat(1:N+1:end) = 0;
    clear Stat
    %---MAT ==========================================
    [~,Stat]=xDF(simts,T,'taper','tukey',sqrt(T),TVFlag);
    Mat.adjt1_amat = Stat.z.rzf;
    Mat.adjt1_amat(1:N+1:end) = 0;
    clear Stat

    [~,Stat]=xDF(simts,T,'taper','tukey',2*sqrt(T),TVFlag);
    Mat.adjt2_amat = Stat.z.rzf;
    Mat.adjt2_amat(1:N+1:end) = 0;
    clear Stat

    %---MAT ==========================================
    [~,Stat]=xDF(simts,T,'taper','curb',T/4,TVFlag);
    Mat.adjc4_amat = Stat.z.rzf;
    Mat.adjc4_amat(1:N+1:end) = 0;
    clear Stat

    %---GLOBALS---------------------------------------------------------------------------------
    disp('Global corrections...')
    %---FSL ==========================================
    [~,adjFSL_amat] = AR1MC(simts,T);
    Mat.adjFSL_amat = adjFSL_amat;
    Mat.adjFSL_amat(1:N+1:end) = 0;
    %---

    %---FOX ==========================================
    Mat.adjFox_amat = FoxBCF(simts,T);
    Mat.adjFox_amat(1:N+1:end) = 0;
    %---

    clear simts
    disp('Generate random networks...')
% RANDOM NETWORKS =============================================================================
    %simts = corrautocorr_sqrtm(zeros(1,N),eye(N),eye(T),T)';
    simts = mvnrnd(zeros(1,N),eye(N),T);
    clear CovMats 

    disp('GM Anal')

    %---MAT ==========================================
    Mat.rmat = atanh(corr(simts)).*sqrt(T-3);
    Mat.rmat(1:N+1:end) = 0;
    %---

    clear simts    
% SAVE ========================================================================================
save(['R/RndNetworks_Yeo_FPP_' num2str(SubID) '_' num2str(j) '_' TVFlag '_Mats.mat'],'Mat')

disp('DONE!')
