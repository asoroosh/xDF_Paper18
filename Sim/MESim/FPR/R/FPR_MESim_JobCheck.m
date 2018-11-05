nRlz=2000;


T_list   = [100, 200, 600, 1200];
rho_list = 0;

for t_cnt=1:numel(T_list)
    for rho_cnt=1:numel(rho_list)
            rt_save=['/storage/essicd/data/HCP/Soroosh/ACAnal/MESim/FPR/R_FPR_MESim_t' num2str(t_cnt) '_r' num2str(rho_cnt) '/'];
            ii = 1; 
	    FJ = [];
            for i_cnt=1:nRlz
                %-----LOAD
                ResStr=[rt_save 'MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(i_cnt) '.mat'];
                
		%exist(ResStr,'file')
		if exist(ResStr,'file')~=2
                    disp(['t' num2str(t_cnt) '_r' num2str(rho_cnt) ' -- ' num2str(i_cnt)])
		    FJ(ii)=i_cnt; 
                    ii=ii+1;	
                end
            end
	save(['FPR_MESim_FailedJobs_t' num2str(t_cnt) '_r' num2str(rho_cnt) '.mat'],'FJ')
	clear FJ
    end
end
