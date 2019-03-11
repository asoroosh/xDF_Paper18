clear
addpath(genpath('~/Home/GitClone/xDF'))

T = 1000; 

nRlz = 1000; 

AC_list  = {0,0.5,[0.8 0.4 0.1],[0.9:-.1:.1]};
rho_list = [0,   0.2, 0.5, 0.7 ,0.9];


% ================ X and Y are indepedent white noises:
CovMat0 = MakeMeCovMat(AC_list{1},T);
CovMat1 = MakeMeCovMat(AC_list{1},T);
CovMat  = cat(3,CovMat0,CovMat1);   
%------------------------------Generate TS

for i = 1:nRlz
    ts_tmp = corrautocorr_sqrtm([0 0],rho_list(1),CovMat,T);
    c_tmp(i)  = corr(ts_tmp(1,:)',ts_tmp(2,:)');
end

histogram(c_tmp,'Normalization','probability')