clear
%close all

%i.e. negative bias means UNDER estimated
%     positive bias means OVER estimated
% E : estimated (somehow)
% C : Calculated (i.e. your ground-truth)
bias = @(MC_C,ParamE) (ParamE-MC_C)./MC_C.*100; 

%---------ANNON FUNCs---------------
%bias = @(old,new) (old-new)./old.*100; %Tom's version
r2z  = @(r,d) atanh(real(r)).*sqrt(d);
r2t  = @(r,d) real(r).*sqrt((d-2)./(1-real(r).^2));
%-----------------------------------

lw=1.2;
fs=15;

nRlz = 5000;

%T_list   = [100, 200, 600, 1200 1500 2000];
T_list   = [100, 200, 600, 1200 1500 2000];
rho_list = [0,   0.2, 0.5, 0.7, 0.9];

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/HetBiv'))

% load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/HetBivStressTest/FC135932/HCP_AC_135932_N47_N87.mat'); % PCC AC profile!

load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/Pocket/RealData_ts_ac.mat','AC_Yeo35','AC_Yeo85','AC_SPPOPGROWDEU')

AC_list  = {0 , 0.5 , AC_SPPOPGROWDEU(2:5) , AC_Yeo35(2:15) , AC_Yeo35(2:20)};%, AC_Yeo85(2:20)};
AC_list_lab={'W','AR1','AR4','AR14','AR20A'};%,'AR20B'};

n_AC_list = numel(AC_list_lab);

for i = 1:n_AC_list
    for ii = 1:n_AC_list
        Labs{i,ii}=[AC_list_lab{i} '-' AC_list_lab{ii}];
    end
end

idx         = find(triu(ones(n_AC_list)));
idx_diag    = find(eye(n_AC_list));
idx_offdiag = setdiff(idx,idx_diag);

%Mrk=[repmat('>',1,11) repmat('o',1,11)]; Mrk=Mrk(randperm(numel(idx_offdiag)));

Mrk=[repmat('>ox',1,15)]; 

spidx={[1 2],[3 4],[5 6],[7 8],[9 10],[11 12]};

fh_bias_SAS  = figure('position',[50,500,750,350]);        hold on; suptitle('Oracle')

%fh_asat_vars = figure('position',[50,500,1800,900]);      hold on; suptitle('Oracle')

fh_bias_CR   = figure('position',[50,500,750,350]);         hold on; suptitle('Oracle')

%fh_z_vars    = figure('position',[50,500,1800,900]);          hold on; suptitle('Variance of Fisher Transformed Corr Coeffs')
%fh_r_vars    = figure('position',[50,500,1800,900]);          hold on; suptitle('Variance of Corr Coeffs')

for t_cnt = 1:numel(T_list)    
    ndp = T_list(t_cnt);
        
    for rho_cnt = 1:numel(rho_list)
        
            %Oracle Simulations------------------------------------------
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
                    disp(['Comb ' num2str(ac_cnt) '  ' AC_list_lab{ac0} ' & ' AC_list_lab{ac1}])
                    %------------------Monster Eq.
                    [ASAt_tmp,CR_tmp] = MonsterEquation_OS(rho_list(rho_cnt),AC_list{ac0},AC_list{ac1},ndp);
                    ASAt(ac0,ac1) = ASAt_tmp;
                    CR(ac0,ac1)   = CR_tmp;
                    clear *_tmp
                    %---------------------------
            end
            %Monte-Carlo Simulations--------------------------------------
            rt_load='/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/R/';
            load([rt_load 'MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(nRlz) '_xDF_TVOff.mat'],'r0')
            
            r0 = r0(1:n_AC_list,1:n_AC_list,:);
            
            %----VARIANCES-----------------------------------------------
            v_r0_tmp = var(r0,[],3);
            v_r0_D  (:,rho_cnt) = v_r0_tmp(idx_diag);
            v_r0_OD (:,rho_cnt) = v_r0_tmp(idx_offdiag);
            
            v_r0_Z_tmp  = var(atanh(r0),[],3);
            v_r0_Z_D (:,rho_cnt)  = v_r0_Z_tmp(idx_diag);
            v_r0_Z_OD (:,rho_cnt) = v_r0_Z_tmp(idx_offdiag);
            
            v_ASAt_diag(:,rho_cnt)    = ASAt(idx_diag);
            v_ASAt_offdiag(:,rho_cnt) = ASAt(idx_offdiag);
            
            %----Biases of V----------------------------------------------
            %--- ASigmaA^\top
            bv_SAS_tmp = real(bias(sqrt(var(r0,[],3)),sqrt(ASAt)));
            
            bv_SAS_diag(:,rho_cnt)    = bv_SAS_tmp(idx_diag);
            bv_SAS_offdiag(:,rho_cnt) = bv_SAS_tmp(idx_offdiag);
            
            all_bv_SAS_diag{t_cnt} = bv_SAS_diag;
            all_bv_SAS_diag{t_cnt} = bv_SAS_offdiag;
            %--- CR
            bv_CR_tmp = bias(sqrt(var(r0,[],3)),sqrt(CR));
            bv_CR_diag(:,rho_cnt)    = bv_CR_tmp(idx_diag);
            bv_CR_offdiag(:,rho_cnt) = bv_CR_tmp(idx_offdiag);
            
            all_bv_CR_diag{t_cnt} = bv_CR_diag;
            all_bv_CR_diag{t_cnt} = bv_CR_offdiag;
            
            clear *_tmp
    end
%     
%     figure(fh_r_vars);
%     sph00=subplot(2,4,spidx{t_cnt}(1)); 
%     hold on; box on; grid on; 
%     title(['T=' num2str(ndp) ', $\lambda_x=\lambda_y$' ],'Interpreter','latex')
%     
%     sph00.XTick      = 1:numel(rho_list);
%     sph00.XTickLabel = {'0','.2','.5','.7','.9'};
%     
%     ylim([0 0.08]); 
%     xlim([1 numel(rho_list)]);
%     
%     for i=1:numel(idx_diag)
%         plot(v_r0_D(i,:),'LineWidth',lw,'marker',Mrk(i));
%     end
%     
%     plot(((1-rho_list.^2).^2)./ndp,'linewidth',lw+.5,'color',[.5 .5 .5],'marker','x') %theoritical variance
%     
%     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
%     ylabel ('$\sigma^2_r$','Interpreter','latex','fontsize',fs)
%     
% %     if t_cnt==4
% %         legend(Labs{idx_diag});    
% %     end
%     %------------------------------------------
%     sph01 = subplot(2,4,spidx{t_cnt}(2)); 
%     hold on; box on; grid on; 
%     title(['T=' num2str(ndp) ', $\lambda_x \neq \lambda_y$' ],'Interpreter','latex')
%     
%     sph01.XTick      = 1:numel(rho_list);
%     sph01.XTickLabel = {'0','.2','.5','.7','.9'};
%     
%     ylim([0 0.08]);
%     xlim([1 numel(rho_list)]);
%     
%     for i=1:numel(idx_offdiag)
%         plot(v_r0_OD(i,:),'LineWidth',lw,'marker',Mrk(i))
%     end
%     
%     plot(((1-rho_list.^2).^2)./ndp,'linewidth',lw+.5,'color',[.5 .5 .5],'marker','x') %theoritical variance
%     
%     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
%     ylabel ('$\sigma^2_r$','Interpreter','latex','fontsize',fs)
    
%     if t_cnt==4
%         legend(Labs{idx_offdiag});
%     end
    %------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------    
%     figure(fh_z_vars);
%     sph00=subplot(2,4,spidx{t_cnt}(1)); 
%     hold on; box on; grid on; 
%     title(['T=' num2str(ndp) ', $\lambda_x=\lambda_y$' ],'Interpreter','latex')
%     
%     sph00.XTick      = 1:numel(rho_list);
%     sph00.XTickLabel = {'0','.2','.5','.7','.9'};
%     
%     ylim([0 0.08]); 
%     xlim([1 numel(rho_list)]);
%     
%     for i=1:numel(idx_diag)
%         plot(v_r0_Z_D(i,:),'LineWidth',lw,'marker',Mrk(i));
%     end
%     
%     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
%     ylabel ('$\sigma^2_z$','Interpreter','latex','fontsize',fs)
%     
% %     if t_cnt==4
% %         legend(Labs{idx_diag});    
% %     end
%     %------------------------------------------    
%     sph01 = subplot(2,4,spidx{t_cnt}(2)); 
%     hold on; box on; grid on; 
%     title(['T=' num2str(ndp) ', $\lambda_x \neq \lambda_y$' ],'Interpreter','latex')
%     
%     sph01.XTick      = 1:numel(rho_list);
%     sph01.XTickLabel = {'0','.2','.5','.7','.9'};
%     
%     ylim([0 0.08]); 
%     xlim([1 numel(rho_list)]);
%     
%     for i=1:numel(idx_offdiag)
%         plot(v_r0_Z_OD(i,:),'LineWidth',lw,'marker',Mrk(i))
%     end
%     
%     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
%     ylabel ('$\sigma^2_z$','Interpreter','latex','fontsize',fs)
    
%     if t_cnt==4
%         legend(Labs{idx_offdiag});
%     end
    %------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------    
    
%     figure(fh_bias_CR);
%     sph00=subplot(2,4,spidx{t_cnt}(1)); 
%     hold on; box on; grid on; 
%     title(['T=' num2str(ndp) ', $\lambda_x=\lambda_y$' ],'Interpreter','latex')
%     
%     sph00.XTick      = 1:numel(rho_list);
%     sph00.XTickLabel = {'0','.2','.5','.7','.9'};
%     
%     ylim([0 0.08]); 
%     xlim([1 numel(rho_list)]);
%     
%     for i=1:numel(idx_diag)
%         plot(v_ASAt_diag(i,:),'LineWidth',lw,'marker',Mrk(i));
%     end
%     
%     plot(((1-rho_list.^2).^2)./ndp,'linewidth',lw+.5,'color',[.5 .5 .5],'marker','x') %theoritical variance
%     
%     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
%     ylabel ('$\sigma^2_r$','Interpreter','latex','fontsize',fs)
%     
% %     if t_cnt==4
% %         legend(Labs{idx_diag});    
% %     end
%     %------------------------------------------    
%     sph01 = subplot(2,4,spidx{t_cnt}(2)); 
%     hold on; box on; grid on; 
%     title(['T=' num2str(ndp) ', $\lambda_x \neq \lambda_y$' ],'Interpreter','latex')
%     
%     sph01.XTick      = 1:numel(rho_list);
%     sph01.XTickLabel = {'0','.2','.5','.7','.9'};
%     
%     ylim([0 0.08]);
%     xlim([1 numel(rho_list)]);
%     
%     for i=1:numel(idx_offdiag)
%         plot(v_ASAt_offdiag(i,:),'LineWidth',lw,'marker',Mrk(i))
%     end
%     
%     plot(((1-rho_list.^2).^2)./ndp,'linewidth',lw+.5,'color',[.5 .5 .5],'marker','x') %theoritical variance
%     
%     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
%     ylabel ('$\sigma^2_r$','Interpreter','latex','fontsize',fs)
%     
%     if t_cnt==4
%         legend(Labs{idx_offdiag});
%     end
    %------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------
    figure(fh_bias_SAS);
    sph10=subplot(2,numel(T_list),spidx{t_cnt}(1)); 
    hold on; box on; grid on; 
    title(['T=' num2str(ndp) ', $\rho_{xx}=\rho_{yy}$' ],'Interpreter','latex')
    
    sph10.XTick      = 1:numel(rho_list);
    sph10.XTickLabel = {'0','.2','.5','.7','.9'};
    sph10.YTick = -50:20:50;
    
    ylim([-50 50]); xlim([1 numel(rho_list)]);
    
    for i=1:numel(idx_diag)
        plot(bv_SAS_diag(i,:),'LineWidth',lw,'marker',Mrk(i));
    end
    
    if t_cnt>3
        xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    end
    %ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    
%     if t_cnt==4
         %legend(Labs{idx_diag},'location','southwest');
%     end
    %------------------------------------------    
    sph11 = subplot(2,numel(T_list),spidx{t_cnt}(2)); 
    hold on; box on; grid on; 
    title(['T=' num2str(ndp) ', $\rho_{xx} \neq \rho_{yy}$' ],'Interpreter','latex')
    
    sph11.XTick      = 1:numel(rho_list);
    sph11.XTickLabel = {'0','.2','.5','.7','.9'};
    sph11.YTick = -50:20:50;
    ylim([-50 50]); xlim([1 numel(rho_list)]);
    for i=1:numel(idx_offdiag)
        plot(bv_SAS_offdiag(i,:),'LineWidth',lw,'marker',Mrk(i))
    end
    
    if t_cnt>3
        xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    end
    %ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    
%     if t_cnt==4
         %legend(Labs{idx_offdiag},'location','southwest');
%     end
    %------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------
    figure(fh_bias_CR);
    sph10=subplot(2,numel(T_list),spidx{t_cnt}(1)); 
    hold on; box on; grid on; 
    title(['T=' num2str(ndp) ', $\rho_{xx}=\rho_{yy}$' ],'Interpreter','latex')
    
    sph10.XTick      = 1:numel(rho_list);
    sph10.XTickLabel = {'0','.2','.5','.7','.9'};
    sph10.YTick = -50:20:50;
    
    ylim([-50 50]); xlim([1 numel(rho_list)]);
    
    for i=1:numel(idx_diag)
        plot(bv_CR_diag(i,:),'LineWidth',lw,'marker',Mrk(i));
    end
    
    if t_cnt>3
        xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    end
    %ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    
%     if t_cnt==4
         %legend(Labs{idx_diag},'location','southwest');
%     end
    %------------------------------------------    
    sph11 = subplot(2,numel(T_list),spidx{t_cnt}(2)); 
    hold on; box on; grid on; 
    title(['T=' num2str(ndp) ', $\rho_{xx} \neq \rho_{yy}$' ],'Interpreter','latex')
    
    sph11.XTick      = 1:numel(rho_list);
    sph11.XTickLabel = {'0','.2','.5','.7','.9'};
    sph11.YTick = -50:20:50;
    ylim([-50 50]); xlim([1 numel(rho_list)]);
    for i=1:numel(idx_offdiag)
        plot(bv_CR_offdiag(i,:),'LineWidth',lw,'marker',Mrk(i))
    end
    
    if t_cnt>3
        xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    end
    %ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    
%     if t_cnt==4
         %legend(Labs{idx_offdiag},'location','southwest');
%     end
    %------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------    
    

end

save('OracleResults.mat','bv_CR_diag','bv_CR_offdiag','bv_SAS_diag','bv_SAS_offdiag');

% set(fh_r_vars,'Color','w'); 
% export_fig(fh_r_vars, 'Figs/Oracle/ME_r_vars.pdf')
% 
% set(fh_z_vars,'Color','w'); 
% export_fig(fh_z_vars, 'Figs/Oracle/ME_z_vars.pdf')
% 
%set(fh_asat_vars,'Color','w'); 
%export_fig(fh_asat_vars,'Figs/Oracle/ME_ASAt.pdf')
% % 
set(fh_bias_CR,'Color','w'); 
%export_fig(fh_bias_CR,'Figs/ME_bias_CnR.pdf')

set(fh_bias_SAS,'Color','w'); 
%export_fig(fh_bias_SAS,'Figs/ME_bias_ASAt.pdf')



%set(fh_bias_wd,'Color','w');  export_fig(fh_bias_wd, 'Figs/ME_bias_bcf_wd.pdf')
%set(fh_bias_wdt,'Color','w'); export_fig(fh_bias_wdt,'Figs/ME_bias_bcf_wdt.pdf')


