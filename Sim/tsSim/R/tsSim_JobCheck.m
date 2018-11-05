nRlz=1000;


T_list   = [100, 200, 600, 1200];
rho_list = [0,   0.2, 0.5, 0.7, 0.9];

for t_cnt=1:4
    for rho_cnt=1:5
            rt_save=['/storage/essicd/data/HCP/Soroosh/ACAnal/tsSim/R_tsSim_t' num2str(t_cnt) '_r' num2str(rho_cnt) '/'];
            for i_cnt=1:nRlz
                %-----LOAD
                ResStr=[rt_save 'tsSim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(i_cnt) '.mat'];
                
%		ResStr

		%exist(ResStr,'file')
		if exist(ResStr,'file')~=2
                    disp(['t' num2str(t_cnt) '_r' num2str(rho_cnt) ' -- ' num2str(i_cnt)])
                end
            end
    end
end
