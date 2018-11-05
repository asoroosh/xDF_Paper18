%PC test:----------------------------------------
%t_cnt=4;
%rho_cnt=4;
%load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/HetBivStressTest/FC135932/HCP_AC_135932_N47_N87.mat','xAC')
%rt_save=[''];
%i_cnt=1;
%------------------------------------------------

addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/DVARS'))

%---------------------------------------------

i_cnt  = str2double(getenv('SGE_TASK_ID'));
rng_sd = randi(i_cnt);
rng(rng_sd)

%---------------------------------------------

% Node 47 (LH-PCC) and 87 (RH-SalVentAttnB_Cinga)
% PCC has a very huge AC, however SalVentAttnB_Cinga is almost zero!
load('/home/wmrnaq/ACAnal/Sim/HetBivStressTest/S/HCP_AC_135932_N47_N87.mat','xAC')
rt_save=['/storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/R_tsSim_t' num2str(t_cnt) '_r' num2str(rho_cnt) '/'];
%---------------------------------------------
T_list   = [100, 200, 600, 1200];
T        = T_list(t_cnt);

rho_list = [0,   0.2, 0.5, 0.7 ,0.9];
rho	 = rho_list(rho_cnt);

AC_list  = {0,0.5,[0.8 0.4 0.1],[0.9:-.1:.1],xAC(1,2:T)};
%---------------------------------------------
    disp('===========')
    disp([num2str(T_list(t_cnt)) '--' num2str(rho_list(rho_cnt))])
    disp('===========')
%---------------------------------------------
nAC     = numel(AC_list);
ncAC    = nAC*(nAC-1)/2+nAC;
[xx,yy] = find(triu(ones(nAC)));
%---------------------------------------------
for ac_cnt=1:ncAC
        ac0 = xx(ac_cnt);
        ac1 = yy(ac_cnt);

        disp(['Comb ' num2str(ac_cnt) '  ' num2str(ac0) ' & ' num2str(ac1)])

        CovMat0 = MakeMeCovMat(AC_list{ac0},T);
        CovMat1 = MakeMeCovMat(AC_list{ac1},T);
        CovMat  = cat(3,CovMat0,CovMat1);   
%------------------------------%------------------------------%------------------------------
        %------------------------------Generate TS
        ts_tmp = corrautocorr_sqrtm([0 0],rho,CovMat,T);
        %-----------------CORR COEFF
        r_tmp = corr(ts_tmp');
        r0_sqrtm(ac0,ac1) = r_tmp(1,2);
        %-----------------Autocorr
        xAC_sqrtm{ac0,ac1} = AC_fft(ts_tmp,T);
        clear *_tmp
%------------------------------%------------------------------%------------------------------
        %------------------------------Generate TS
        ts_tmp = corrautocorr([0 0],rho,CovMat,T);
        clear CovMat
        %-----------------CORR COEFF
        r_tmp = corr(ts_tmp');
        r0_chol(ac0,ac1) = r_tmp(1,2);
        %-----------------Autocorr
        xAC_chol{ac0,ac1} = AC_fft(ts_tmp,T);
        clear *_tmp CovMat

end



    if exist(rt_save,'dir')~=7; mkdir(rt_save); end
    
    save([rt_save 'tsSim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(i_cnt) '.mat'],...
        'r0_sqrtm',...
        'r0_chol',...
        'xAC_sqrtm',...
        'xAC_chol',...
        'rng_sd');

    disp('DONE!')


