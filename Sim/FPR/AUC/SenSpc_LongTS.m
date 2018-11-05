addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/DVARS'))
addpath(genpath('~/bin/2017_01_15_BCT'))

i_cnt  = str2double(getenv('SGE_TASK_ID'));

%------------------------------INITIALISE--------------------------------------------
rho_list = [0 0.2 0.5 0.7 0.9];
T_list   = [1500 2000]; %Oi! before adding 1500 and 2000, becareful that the size(mts,1) is only 1200

Tp = T_list(t_cnt);

TarDen = 0.15; %target density

TVFlag = 'TVOn';

rho_u = 9.5;
rho_d = 8.5;

upb   = rho_u;
downb = rho_d;

Np = 114;

rng_sd = randi(i_cnt);
rng(rng_sd)

%load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/S/HCP_MPP_135932_ROIsTS.mat'],'mts');
%load(['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FC/Yeo/FPP/HCP_FPP_135932_ROIsTS.mat'],'mts')

SubID = 118528; %This was found based on the autocorrelation lengths calculated via global BCF
%SubID = 135932;
load(['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FC/Yeo/FPP/HCP_FPP_' num2str(SubID) '_ROIsTS.mat'],'mts')

%T  = size(mts,1);
%N  = size(mts,2);

disp(['T: ' num2str(Tp) ])
disp(['RhoU: ' num2str(upb)  ' _ RhoD: ' num2str(downb) ])

%--------------------------------BODY-------------------------------------------------
if Tp>1200
	Tpp = 1200;
end

AC = AC_fft(mts(1:Tpp,:),Tpp); 
AC = AC(:,2:end-1);

aCovMat = zeros(Tp,Tp,Np);
for i=1:Np
	i
	aCovMat(:,:,i) = MakeMeCovMat(AC(i,:),Tp);
end

%----This is for a synthetic adj matrix------------------------------------------------
%blkrndrng = @(downb,upb,NN) ((upb-downb).*rand(NN,NN) + downb)./10;
%rCovMat = blkdiag(blkrndrng(downb,upb,blksz),blkrndrng(downb,upb,blksz),...
%    blkrndrng(downb,upb,blksz),blkrndrng(downb,upb,blksz),blkrndrng(downb,upb,blksz)...
%    ,blkrndrng(downb,upb,blksz),blkrndrng(downb,upb,blksz),blkrndrng(downb,upb,(Np-blksz*7)));

% WhatDiags=[50 90];
% for didx=WhatDiags
%     rCovMat = rCovMat + diag(vecrndrng(downb,upb,didx),Np-didx) + diag(vecrndrng(downb,upb,didx),Np-didx)';
% end
%rCovMat(30:30+blksz-1,70:70+blksz-1) = blkrndrng(downb,upb,blksz);

%rCovMat = rCovMat+fliplr(rCovMat);
%rCovMat(16*3+1:16*3+16,16*3+1:16*3+16) = rCovMat(16*3+1:16*3+16,16*3+1:16*3+16)./2; 

%----This is for real data------------------------------------------------------------
mat = corr(mts);
mat (mat<0) = 0;
mat = threshold_proportional(mat,TarDen);
mat (mat>0) =1;

Idx = find(mat(1:Np,1:Np));
%----------------------------------CovMat--------------------------------------------
vecrndrng = @(downb,upb,NN) ((upb-downb).*rand(NN,1) + downb)./10;

rCovMat = zeros(Np,Np);
rCovMat(Idx) = vecrndrng(downb,upb,numel(Idx));

rCovMat(1:Np+1:end) = 1; %temprorily set it to one to generate the time series. 
disp('covmats generated.')

%----------------------------------Simulate TS---------------------------------------
simts = corrautocorr(zeros(1,Np),rCovMat,aCovMat,Tp);
disp('ts generated.')
size(simts)
%-------------------------FIND POS & NEGS -------------------------------------------
rCovMat(1:Np+1:end) = 0; %set back the diagonal to zero, because we don't care about it in Sen Spc analysis.
CG = Np*(Np-1)./2;
POS = numel(find(rCovMat));
NEG = CG-POS;
%-----------------------------------------------------------------
%----ME
[V,Stat]      = xDF(simts,Tp,'taper','shrink',TVFlag);
p_MEs1 = Stat.p.f_Pval;
Z_MEs1 = Stat.z.rzf;
clear Stat 
%----ME
[V,Stat]      = xDF(simts,Tp,'taper','tukey',2*sqrt(Tp),TVFlag);
p_MEt2 = Stat.p.f_Pval;
Z_MEt2 = Stat.z.rzf;
clear Stat 
%----ME
[V,Stat]      = xDF(simts,Tp,TVFlag);
p_ME = Stat.p.f_Pval;
Z_ME = Stat.z.rzf;
clear Stat 
%----ME
[V,Stat]      = xDF(simts,Tp,'taper','tukey',sqrt(Tp),TVFlag);
p_MEt1 = Stat.p.f_Pval;
Z_MEt1 = Stat.z.rzf;
clear Stat 
%----ME
[V,Stat]      = xDF(simts,Tp,'taper','curb',Tp/4,TVFlag);
p_MEc4 = Stat.p.f_Pval;
Z_MEc4 = Stat.z.rzf;
clear Stat 
%------------------------------
%----CR
[~,Z_CR,p_CR]       = CRBCF(simts,Tp);
%----CH
[~,Z_CH,p_CH]       = CheltonBCF(simts,Tp);
%----AR1
[~,Z_AR1,p_AR1]     = Bartlett46_fft(simts,Tp);
%----AR1MC
[~,Z_AR1MC,p_AR1MC] = AR1MC(simts,Tp);
%----Fox
[Z_Fox,~,p_Fox]     = FoxBCF(simts,Tp);
%----Naive
[R_Naive,p_Naive]   = corr(simts');
Z_Naive = atanh(R_Naive).*sqrt(Tp-3);

EstsLables = {'ME','MEs1','MEt1','MEt2','MEc4','CR','CH','AR1','AR1MC','Fox','Naive'};

for m_cnt=1:numel(EstsLables)
    disp(['Z_' EstsLables{m_cnt}])
    
    %P_tmp = eval(['p_' EstsLables{m_cnt}]);
    %H_tmp = fdr_bh(P_tmp);

   %------------------------------------------------------------------    
    %TP_tmp = numel(find(H_tmp.*rCovMat));
    %FP_tmp = numel(find(H_tmp.*~rCovMat));
    %FN_tmp = numel(find(~H_tmp.*rCovMat));
    %TN_tmp = NEG-FP_tmp;

    %Sen(1,m_cnt) = TP_tmp./(TP_tmp+FN_tmp);
    %Spc(1,m_cnt) = TN_tmp./(TN_tmp+FP_tmp);
    %Acc(1,m_cnt) = (TP_tmp+TN_tmp)./(TP_tmp+FP_tmp+FN_tmp+TN_tmp);
    %------------------------------------------------------------------ 
    
    %P(:,:,m_cnt) = P_tmp;
    %H(:,:,m_cnt) = H_tmp;
    
    Z(:,:,m_cnt) = eval(['Z_' EstsLables{m_cnt}]);

    %TP(1,m_cnt) = TP_tmp;
    %FP(1,m_cnt) = FP_tmp;
    %FN(1,m_cnt) = FN_tmp;
    %TN(1,m_cnt) = TN_tmp;
    
    clear *_tmp
end

%----SAVE
save(['R/Sen_Spec_t' num2str(Tp) '_' num2str(SubID) '_' num2str(i_cnt) '_' TVFlag '_r' num2str(rho_u) '-' num2str(rho_d) '_JustZs.mat'],'Z','rCovMat');
disp('Done!')
%---DONE

