%PC test:----------------------------------------
%t_cnt=4;
%rho_cnt=4;
%load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/HetBivStressTest/FC135932/HCP_AC_135932_N47_N87.mat','xAC')
%rt_save=[''];
%i_cnt=1;
%------------------------------------------------

addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/DVARS'))

i_cnt  = str2double(getenv('SGE_TASK_ID'));

TVFlag = 'TVOff';

%------------------------------------THIS SHOULD BE REMOVED!-------------------------------
%load(['/home/wmrnaq/ACAnal/Sim/MESim/R/MESim_FailedJobs_t' num2str(t_cnt) '_r' num2str(rho_cnt) '.mat'],'FJ')

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

rt_save=['/storage/essicd/data/HCP/Soroosh/ACAnal/MESim/R_MESim_t' num2str(t_cnt) '_r' num2str(rho_cnt) '/'];

T_list   = [100, 200, 600, 1200, 1500, 2000];
T        = T_list(t_cnt);

rho_list = [0,   0.2, 0.5, 0.7 ,0.9];
rho      = rho_list(rho_cnt);

load('/home/wmrnaq/ACAnal/Sim/MESim/S/RealData_ts_ac.mat','AC_Yeo35','AC_Yeo85','AC_SPPOPGROWDEU')

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
        %------------------Monster Eq.----------------------------------------
        [ASAt_tmp,Stat_tmp] = xDF(ts,T,TVFlag);
        
	ASAt(ac0,ac1) = ASAt_tmp(1,2);

        p_f(ac0,ac1)  = Stat_tmp.p.f_Pval(1,2);
	z_f(ac0,ac1)  = Stat_tmp.z.rzf(1,2);
%	p_r(ac0,ac1)  = Stat_tmp.p.r_Pval(1,2);

	%------------------Clifford & Richardson w/out Curbing
        
%	CR(ac0,ac1)   = Stat_tmp.CR;	

	clear *_tmp
        %---------------------------
        %------------------Monster Eq. Tukey Tapering -- Woolrich
        [ASAt_tmp,Stat_tmp] = xDF(ts,T,'taper','tukey',round(sqrt(T)),TVFlag);

        ASAt_tt1(ac0,ac1) = ASAt_tmp(1,2);

        p_f_tt1(ac0,ac1)  = Stat_tmp.p.f_Pval(1,2);
 %       p_r_tt1(ac0,ac1)  = Stat_tmp.p.r_Pval(1,2);

	z_f_tt1(ac0,ac1)  = Stat_tmp.z.rzf(1,2);
 %       z_r_tt1(ac0,ac1)  = Stat_tmp.z.rz(1,2);

        clear *_tmp
        %---------------------------
        %------------------Monster Eq. Tukey Tapering -- Chatfield
        [ASAt_tmp,Stat_tmp] = xDF(ts,T,'taper','tukey',2*round(sqrt(T)),TVFlag);

        ASAt_tt2(ac0,ac1) = ASAt_tmp(1,2);

        p_f_tt2(ac0,ac1)  = Stat_tmp.p.f_Pval(1,2);
 %       p_r_tt2(ac0,ac1)  = Stat_tmp.p.r_Pval(1,2);

	z_f_tt2(ac0,ac1)  = Stat_tmp.z.rzf(1,2);
 %       z_r_tt2(ac0,ac1)  = Stat_tmp.z.rz(1,2);

        clear *_tmp
        %---------------------------
        %------------------Monster Eq. Curb Tapering -- Anderson
        [ASAt_tmp,Stat_tmp] = xDF(ts,T,'taper','curb',round(T/4),TVFlag);

        ASAt_tc4(ac0,ac1) = ASAt_tmp(1,2);

        p_f_tc4(ac0,ac1)  = Stat_tmp.p.f_Pval(1,2);
%        p_r_tc4(ac0,ac1)  = Stat_tmp.p.r_Pval(1,2);

	z_f_tc4(ac0,ac1)  = Stat_tmp.z.rzf(1,2);
%        z_r_tc4(ac0,ac1)  = Stat_tmp.z.rz(1,2);

        clear *_tmp
        %---------------------------
        %------------------Monster Eq. Curb Tapering -- Pyper 
        [ASAt_tmp,Stat_tmp] = xDF(ts,T,'taper','curb',round(T/5),TVFlag);

        ASAt_tc5(ac0,ac1) = ASAt_tmp(1,2);

        p_f_tc5(ac0,ac1)  = Stat_tmp.p.f_Pval(1,2);
%        p_r_tc5(ac0,ac1)  = Stat_tmp.p.r_Pval(1,2);

        z_f_tc5(ac0,ac1)  = Stat_tmp.z.rzf(1,2);
%        z_r_tc5(ac0,ac1)  = Stat_tmp.z.rz(1,2);

        clear *_tmp
        %---------------------------
        %------------------Monster Eq. Shrink Tapering 1 period-- Afyouni & Nichols 
        [ASAt_tmp,Stat_tmp] = xDF(ts,T,'taper','shrink',1,TVFlag);

        ASAt_ts1(ac0,ac1) = ASAt_tmp(1,2);

        p_f_ts1(ac0,ac1)  = Stat_tmp.p.f_Pval(1,2);
%        p_r_ts1(ac0,ac1)  = Stat_tmp.p.r_Pval(1,2);

        z_f_ts1(ac0,ac1)  = Stat_tmp.z.rzf(1,2);
%        z_r_ts1(ac0,ac1)  = Stat_tmp.z.rz(1,2);

        clear *_tmp
        
	%---------------------------
	%-- Chelthon with curbing factor T/5; exactly what Pyper suggests
	[v_ch_tmp,z_ch_tmp,p_ch_tmp] = CheltonBCF(ts,T);

	v_ch(ac0,ac1) = v_ch_tmp(1,2);
	z_ch(ac0,ac1) = z_ch_tmp(1,2);
	p_ch(ac0,ac1) = p_ch_tmp(1,2);

	%---------------------------
        %-- CR with curbing factor T/5; exactly what Pyper suggests
        [v_cr_tmp,z_cr_tmp,p_cr_tmp] = CRBCF(ts,T);

        v_cr(ac0,ac1) = v_cr_tmp(1,2);
        z_cr(ac0,ac1) = z_cr_tmp(1,2);
        p_cr(ac0,ac1) = p_cr_tmp(1,2);


	%---------------------------
        %-- CR with curbing factor T/5; exactly what Pyper suggests
        [v_ar1_tmp,z_ar1_tmp,p_ar1_tmp] = Bartlett46_fft(ts,T);

        v_ar1(ac0,ac1) = v_ar1_tmp(1,2);
        z_ar1(ac0,ac1) = z_ar1_tmp(1,2);
        p_ar1(ac0,ac1) = p_ar1_tmp(1,2);
	
	clear *_tmp
	
end
    
    if exist(rt_save,'dir')~=7; mkdir(rt_save); end
    
    save([rt_save 'MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(i_cnt) '_xDF_' TVFlag '.mat'],...
        'r0','p_naive',...
        'ASAt','p_f','z_f',...
	'ASAt_tt1','p_f_tt1','z_f_tt1',...
	'ASAt_tt2','p_f_tt2','z_f_tt2',...
	'ASAt_tc4','p_f_tc4','z_f_tc4',... 
	'ASAt_tc5','p_f_tc5','z_f_tc5',...
	'ASAt_ts1','p_f_ts1','z_f_ts1',... %'ASAt_ts2','p_f_ts2','p_r_ts2','z_f_ts2','z_r_ts2',...
	'v_ch' ,'z_ch' ,'p_ch',...
	'v_cr' ,'z_cr' ,'p_cr',...
	'v_ar1','z_ar1','p_ar1',...        
	'rng_sd');

    disp('DONE!')

