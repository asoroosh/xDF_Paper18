clear

nRlz=1000;

T_list   = [100, 200, 600, 1200];
rho_list = [0,   0.2, 0.5, 0.7 0.9];

for t_cnt=1:4
    for rho_cnt=1:5
	    disp(['t' num2str(t_cnt) '_r' num2str(rho_cnt)])

            rt_save=['/storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/R_tsSim_t' num2str(t_cnt) '_r' num2str(rho_cnt) '/'];
            for i_cnt=1:nRlz
                
		%-----LOAD
                Res=load([rt_save 'tsSim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(i_cnt) '.mat']);
                
		%-----CORR COEFF
                r0_sqrtm(:,:,i_cnt) = Res.r0_sqrtm;
                r0_chol (:,:,i_cnt) = Res.r0_chol; 
		%-----BCF
                xAC_chol  {i_cnt} = Res.xAC_chol;
                xAC_sqrtm {i_cnt} = Res.xAC_sqrtm;

                clear Res
            end
            save(['R/tsSim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '.mat'],...
                'T_list','rho_list','nRlz',...
                'r0_sqrtm','r0_chol',...
		'xAC_sqrtm','xAC_chol');
            
            clear r0 CF_ON V_ON EDF_ON CF_OFF V_OFF EDF_OFF CF_AR1
    end
end
