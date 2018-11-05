clear

nRlz=2000;

T_list   = [100, 200, 600, 1200];
rho_list = 0;

for t_cnt=1:numel(T_list)
    for rho_cnt=1:numel(rho_list)
	    disp(['t' num2str(t_cnt) '_r' num2str(rho_cnt)])

            rt_save=['/storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/R_FPR_MESim_t' num2str(t_cnt) '_r' num2str(rho_cnt) '/'];
            for i_cnt=1:nRlz
                
		%-----LOAD---------------------------------------
                Res=load([rt_save 'MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(i_cnt) '.mat']);
                
		%-----CORR COEFF---------------------------------
                r0(:,:,i_cnt)   = Res.r0;
                
		%- pval Estimates----------------------------------
		% PARAM
		p_naive    (:,:,i_cnt) = Res.p_naive;
		%p_f       (:,:,i_cnt) = Res.p_f;		
		%p_f_tc4   (:,:,i_cnt) = Res.p_f_tc4;
		%p_f_tc5   (:,:,i_cnt) = Res.p_f_tc5;

		p_f_ts1   (:,:,i_cnt) = Res.p_f_ts1;
		%p_f_ts2   (:,:,i_cnt) = Res.p_f_ts2;

		%p_f_tt1   (:,:,i_cnt) = Res.p_f_tt1;
		%p_f_tt2   (:,:,i_cnt) = Res.p_f_tt2;

		%NON PARAM----------------------------------------		
		%monster eq
		p_np_f_ts1(:,:,i_cnt)  = Res.p_np_f_ts1;
		%fisher
		p_np_fisher(:,:,i_cnt) = Res.p_np_fisher;		

                clear Res
            end
            save(['R/FPR_MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '.mat'],...
                'T_list','rho_list','nRlz',...
                'r0','p_naive',...
		'p_f_ts1','p_np_f_ts1','p_np_fisher');
		%'p_f' ,'p_f_tc4' ,'p_f_tc5' ,'p_f_ts1' ,'p_f_ts2' ,'p_f_tt1' ,'p_f_tt2'

            clear r0 CR ASAt ASAt_* z_* p_* Res
    end
end
