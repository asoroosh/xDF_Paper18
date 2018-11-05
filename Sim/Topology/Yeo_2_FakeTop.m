clear

j = str2double(getenv('SGE_TASK_ID'));

rng(j)

thrng=0.01:0.01:0.50;

%SubID = 135932;
SubID = 118528;

TVFlag = 'TVOff';

addpath(genpath('~/bin/2017_01_15_BCT'))
addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/contest'))

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
    amat = atanh(corr(simts')).*sqrt(T-3);
    amat(1:N+1:end) = 0;
    %---
% ADJUSTED NETWORKS ==================================
	disp('AC corrections...')
    %---MAT ==========================================
    [~,ZCR] = CRBCF(simts,T);
    adjCR_amat = ZCR;
    adjCR_amat(1:N+1:end) = 0;
    
    clear ZCR
    %---

    %---MAT ==========================================
    [~,Stat]=xDF(simts,T,TVFlag);
    adj_amat = Stat.z.rzf; 
    adj_amat(1:N+1:end) = 0;
    
    clear Stat
    %---

    %---MAT ==========================================
    [~,Stat]=xDF(simts,T,'taper','shrink',TVFlag);
    adjs1_amat = Stat.z.rzf;
    adjs1_amat(1:N+1:end) = 0;

    clear Stat
    %---

    %---FSL ==========================================
    [~,ZFSL] = AR1MC(simts,T);
    adjFSL_amat = ZFSL;
    adjFSL_amat(1:N+1:end) = 0;
    %---

    %---FOX ==========================================
    ZFox = FoxBCF(simts,T);
    adjFox_amat = ZFox;
    adjFox_amat(1:N+1:end) = 0;
    %---

    clear simts
	disp('Generate random networks...')
% RANDOM NETWORKS =============================================================================

    %simts = corrautocorr_sqrtm(zeros(1,N),eye(N),eye(T),T)';
    simts = mvnrnd(zeros(1,N),eye(N),T);
    clear CovMats 

    disp('GM Anal')

    %---MAT ==========================================
    rmat = atanh(corr(simts)).*sqrt(T-3);
    rmat(1:N+1:end) = 0;
    %---

    clear simts    

disp('thresholding')
% GM ========================================================================================
    thr_cnt=1;
    for thr=thrng
	disp(['thr: ' num2str(thr)])

        %AC NETWORKS ------------------------------------------------------
        t_m_cm_ac = threshold_proportional(abs(amat),thr);
        t_m_cm_ac(t_m_cm_ac>0) = 1; %binarise 

        ACGM.m_le_ac(thr_cnt) = mean(efficiency_bin(t_m_cm_ac,1));
        ACGM.m_ge_ac(thr_cnt) = efficiency_bin(t_m_cm_ac);

	clear t_cm_*
        %AC NETWORKS ------------------------------------------------------
        
	%Without shrinking-------------------------------------------------
	t_m_cm_adjac = threshold_proportional(abs(adj_amat),thr);
        t_m_cm_adjac(t_m_cm_adjac>0) = 1; %binarise 

        ACGM.m_le_adjac(thr_cnt) = mean(efficiency_bin(t_m_cm_adjac,1));
        ACGM.m_ge_adjac(thr_cnt) = efficiency_bin(t_m_cm_adjac);
        
	clear t_cm_*
	%By Shrinkage------------------------------------------------------
	t_m_cm_adjac_s1 = threshold_proportional(abs(adjs1_amat),thr);
        t_m_cm_adjac_s1(t_m_cm_adjac_s1>0) = 1; %binarise 

        ACGM.m_le_adjac_s1(thr_cnt) = mean(efficiency_bin(t_m_cm_adjac_s1,1));
        ACGM.m_ge_adjac_s1(thr_cnt) = efficiency_bin(t_m_cm_adjac_s1);

	clear t_cm_*
	%By Shrinkage------------------------------------------------------
        t_m_cm_adjCR_amat = threshold_proportional(abs(adjCR_amat),thr);
        t_m_cm_adjCR_amat(t_m_cm_adjCR_amat>0) = 1; %binarise 

        ACGM.m_le_adjac_cr(thr_cnt) = mean(efficiency_bin(t_m_cm_adjCR_amat,1));
        ACGM.m_ge_adjac_cr(thr_cnt) = efficiency_bin(t_m_cm_adjCR_amat);

	clear t_cm_*
	%By FSLNets--------------------------------------------------------
        t_m_cm_adjac_fsl = threshold_proportional(abs(adjFSL_amat),thr);
        t_m_cm_adjac_fsl(t_m_cm_adjac_fsl>0) = 1; %binarise 

        ACGM.m_le_adjac_fsl(thr_cnt) = mean(efficiency_bin(t_m_cm_adjac_fsl,1));
        ACGM.m_ge_adjac_fsl(thr_cnt) = efficiency_bin(t_m_cm_adjac_fsl);

	clear t_cm_*
	%By Fox et al 2005-------------------------------------------------
        t_m_cm_adjac_fox = threshold_proportional(abs(adjFox_amat),thr);
        t_m_cm_adjac_fox(t_m_cm_adjac_fox>0) = 1; %binarise 

        ACGM.m_le_adjac_fox(thr_cnt) = mean(efficiency_bin(t_m_cm_adjac_fox,1));
        ACGM.m_ge_adjac_fox(thr_cnt) = efficiency_bin(t_m_cm_adjac_fox);
        
	clear t_cm_*
	%Erdos-Renyi Network ----------------------------------------------
        erm = full(erdrey(N,round(thr*(N*(N-1)/2))));

        ACGM.m_le_erm(thr_cnt) = mean(efficiency_bin(erm,1));
        ACGM.m_ge_erm(thr_cnt) = efficiency_bin(erm);
        
        %RND NETWORKS ------------------------------------------------------
        t_m_cm_rnd = threshold_proportional(abs(rmat),thr);
        t_m_cm_rnd(t_m_cm_rnd>0) = 1;

        ACGM.m_le_rnd(thr_cnt) = mean(efficiency_bin(t_m_cm_rnd,1));
        ACGM.m_ge_rnd(thr_cnt) = efficiency_bin(t_m_cm_rnd);
        
	clear t_cm_*
        %LATTICE ---------------------------------------------------------
        lat = makelatticeCIJ(N,round(thr*(N*(N-1))));

        ACGM.m_le_lat(thr_cnt) = mean(efficiency_bin(lat,1));
        ACGM.m_ge_lat(thr_cnt) = efficiency_bin(lat);
        
        thr_cnt=thr_cnt+1;
    end

% SAVE ========================================================================================
save(['R/RndNetworks_Yeo_FPP_' num2str(SubID) '_' num2str(j) '_' TVFlag '_Top.mat'],'rmat','amat','adj_amat','adjs1_amat','ACGM','adjCR_amat','adjFSL_amat','adjFox_amat')

disp('DONE!')
