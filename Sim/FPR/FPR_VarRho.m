addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/DVARS'))
addpath(genpath('~/bin/2017_01_15_BCT'))

i_cnt  = str2double(getenv('SGE_TASK_ID'));

rho_list = [0 0.2 0.5 0.7 0.9];

rho = rho_list(rho_cnt)*10;

rho_u = rho;
rho_d = rho;

rng_sd = randi(i_cnt);
rng(rng_sd)

%load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/S/HCP_MPP_135932_ROIsTS.mat'],'mts');
load(['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FC/Yeo/FPP/HCP_FPP_135932_ROIsTS.mat'],'mts')

T  = size(mts,1);
N  = size(mts,2);

AC = AC_fft(mts,T); 
AC = AC(:,2:101);

fs = 12; 

Tp = T;

Np = 112;

CG = Np*(Np-1)./2;

aCovMat = zeros(Tp,Tp,Np);
for i=1:Np
    aCovMat(:,:,i) = MakeMeCovMat_sqrtm(AC(i,:),Tp);
end

blkrndrng = @(downb,upb,NN) ((upb-downb).*rand(NN,NN) + downb)./10;
vecrndrng = @(downb,upb,NN) ((upb-downb).*rand(NN,1) + downb)./10;

blksz = 16;
upb   = rho_u;
downb = rho_d;


disp(['RhoU: ' num2str(upb)  ' _ RhoD: ' num2str(downb) ])

%----This is for a synthetic adj matrix
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

%----This is for real data
mat = corr(mts);
mat (mat<0) = 0;
mat = threshold_proportional(mat,.15);
mat (mat>0) =1;

Idx = find(mat(1:Np,1:Np));

rCovMat = zeros(Np,Np);
rCovMat(Idx) = vecrndrng(downb,upb,numel(Idx));

rCovMat(1:Np+1:end) = 1; %temprorily set it to one to generate the time series. 

disp('covmats generated.')

POS = numel(find(rCovMat))-Np;
NEG = CG-POS;

simts = corrautocorr(zeros(1,Np),rCovMat,aCovMat,Tp);
disp('ts generated.')

rCovMat(1:Np+1:end) = 0; %set back the diagonal to zero, because we don't care about it in Sen Spc analysis. 

%----ME
[V,Stat]      = PearCorrVarEst(simts,Tp,'taper','shrink',1);
p_ME = Stat.p.f_Pval;
%----CR
[~,~,p_CR]    = CRBCF(simts,Tp);
%----CH
[~,~,p_CH]    = CheltonBCF(simts,Tp);
%----AR1
[~,~,p_AR1]   = Bartlett46_fft(simts,Tp);
%----AR1MC
[~,~,p_AR1MC] = AR1MC(simts,Tp);
%----Fox
[~,~,p_Fox]   = FoxBCF(simts,Tp);
%----Naive
[~,p_Naive]   = corr(simts');

EstsLables = {'ME','CR','CH','AR1','AR1MC','Fox','Naive'};

for m_cnt=1:numel(EstsLables)
    disp(['p_' EstsLables{m_cnt}])
    
    P_tmp = eval(['p_' EstsLables{m_cnt}]);
    H_tmp = fdr_bh(P_tmp);
    
    if strcmp(EstsLables{m_cnt},'ME')
        H_tmp(1:Np+1:end) = 0;
    end

    if strcmp(EstsLables{m_cnt},'CR')
        H_tmp(1:Np+1:end) = 1;
    end

    TP_tmp = numel(find(H_tmp.*rCovMat));
    FP_tmp = numel(find(H_tmp.*~rCovMat));
    FN_tmp = numel(find(~H_tmp.*rCovMat));
    TN_tmp = NEG-FP_tmp;

    Sen(1,m_cnt) = TP_tmp./(TP_tmp+FN_tmp);
    Spc(1,m_cnt) = TN_tmp./(TN_tmp+FP_tmp);
    Acc(1,m_cnt) = (TP_tmp+TN_tmp)./(TP_tmp+FP_tmp+FN_tmp+TN_tmp);
        
    P(:,:,m_cnt) = P_tmp;
    H(:,:,m_cnt) = H_tmp;

    TP(1,m_cnt) = TP_tmp;
    FP(1,m_cnt) = FP_tmp;
    FN(1,m_cnt) = FN_tmp;
    TN(1,m_cnt) = TN_tmp;
    
    clear *_tmp
end

%----SAVE
save(['R/Sen_Spec_' num2str(rho_cnt) '_' num2str(i_cnt) '.mat'],'Sen','Spc','Acc','H','TP','FP','FN','TN');

disp('Done!')
%---DONE

