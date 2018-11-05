addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/DVARS'))

i_cnt  = str2double(getenv('SGE_TASK_ID'));

%------------------------------------THIS SHOULD BE REMOVED!-------------------------------
%load(['/home/wmrnaq/ACAnal/Sim/MESim/FPR/R/MESim_FailedJobs_t' num2str(t_cnt) '_r' num2str(rho_cnt) '.mat'],'FJ')

%if isempty(FJ) 
%	exit; 
%end

%i_cnt = FJ(i_cnt);
%disp(['This is the failed job number: ' num2str(i_cnt) ' re-running....'])
%------------------------------------THIS SHOULD BE REMOVED!-------------------------------

rng_sd = randi(i_cnt);
rng(rng_sd)

% Node 47 (LH-PCC) and 87 (RH-SalVentAttnB_Cinga)
% PCC has a very huge AC, however SalVentAttnB_Cinga is almost zero!

rt_save=['/storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/R_FPR_MESim_t' num2str(t_cnt) '_r' num2str(rho_cnt) '/'];

T_list   = [100, 200, 600, 1200];
T        = T_list(t_cnt);

rho_list = [0,   0.2, 0.5, 0.7 ,0.9];
rho      = rho_list(rho_cnt);

load('/home/wmrnaq/ACAnal/Sim/MESim/S/RealData_ts_ac.mat','AC_Yeo35','AC_Yeo85','AC_SPPOPGROWDEU')

nRlz = 2000;

AC_list  = {0 , 0.5 , AC_SPPOPGROWDEU(2:5) , AC_Yeo35(2:15) , AC_Yeo35(2:20) , AC_Yeo85(2:20)};

    disp('===========')
    disp([num2str(T_list(t_cnt)) '--' num2str(rho_list(rho_cnt))])
    disp('===========')

nAC     = numel(AC_list);
ncAC    = nAC*(nAC-1)/2+nAC;
[xx,yy] = find(triu(ones(nAC)));
%---------------------------------------------
for ac_cnt=1:ncAC
        ac0 = xx(ac_cnt);
        ac1 = yy(ac_cnt);
        
        disp(['Comb ' num2str(ac_cnt) '  ' num2str(ac0) ' & ' num2str(ac1)])
        
        CovMat0=MakeMeCovMat(AC_list{ac0},T);
        CovMat1=MakeMeCovMat(AC_list{ac1},T);
        CovMat=cat(3,CovMat0,CovMat1);
        %------------------------------Generate TS
        ts=corrautocorr_sqrtm([0 0],rho,CovMat,T);
        clear CovMat
        %-----------------CORR COEFF
        [r_tmp,p_tmp]    = corr(ts');
        r0(ac0,ac1) 	 = r_tmp(1,2);
	p_naive(ac0,ac1) = p_tmp(1,2);
        clear *_tmp
        
	%------------------Monster Eq.
        % [p_np_f(ac0,ac1)     p_f(ac0,ac1)]     = NonParamCorrInf(ts,T,nRlz); 
        %---------------------------
        %------------------Monster Eq. Tukey Tapering -- Woolrich
        % [p_np_f_tt1(ac0,ac1) p_f_tt1(ac0,ac1)] = NonParamCorrInf(ts,T,nRlz,'taper','tukey',round(sqrt(T)));
        %---------------------------
        %------------------Monster Eq. Tukey Tapering -- Chatfield
        % [p_np_f_tt2(ac0,ac1),p_np_fisher(ac0,ac1),p_f_tt2(ac0,ac1)] = NonParamCorrInf(ts,T,nRlz,'taper','tukey',2*round(sqrt(T)));
        %---------------------------
        %------------------Monster Eq. Curb Tapering -- Anderson
        % [p_np_f_tc4(ac0,ac1) p_f_tc4(ac0,ac1)] = NonParamCorrInf(ts,T,nRlz,'taper','curb',round(T/4));
        %---------------------------
        %------------------Monster Eq. Curb Tapering -- Pyper 
        % [p_np_f_tc5(ac0,ac1) p_f_tc5(ac0,ac1)] = NonParamCorrInf(ts,T,nRlz,'taper','curb',round(T/5)); 
        %---------------------------
        %------------------Monster Eq. Shrink Tapering 1 period-- Afyouni & Nichols 
        [p_np_f_ts1(ac0,ac1) p_np_fisher(ac0,ac1) p_f_ts1(ac0,ac1)] = NonParamCorrInf(ts,T,nRlz,'taper','shrink',1); 
        %---------------------------
        %------------------Monster Eq. Shrink Tapering 2 period-- Afyouni & Nichols 
        % p_f_ts2(ac0,ac1)  = NonParamCorrInf(ts,T,nRlz,'taper','shrink',2);
        %---------------------------

end
    
%'p_np_f','p_np_f_tt1','p_np_f_tt2','p_np_f_tc4','p_np_f_tc5',
%'p_f','p_f_tt1','p_f_tt2','p_f_tc4','p_f_tc5'
    if exist(rt_save,'dir')~=7; mkdir(rt_save); end
    
    save([rt_save 'MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(i_cnt) '.mat'],...
        'r0','p_naive',...
	'p_np_f_ts1','p_np_fisher',...
	'p_f_ts1',...
	'rng_sd');

    disp('DONE!')
