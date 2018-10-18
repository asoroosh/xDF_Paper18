clear
%close all

%---------ANNON FUNCs---------------
%i.e. negative bias means UNDER estimated
%     positive bias means OVER estimated
% E : estimated (i.e. Monster Equation)
% C : Calculated (i.e. your ground-truth; Monte Carlo simulations)
bias = @(MC_C,ParamE) (ParamE-MC_C)./MC_C.*100;

%bias = @(old,new) (new-old)./old.*100; %Tom's version!
r2z  = @(r,d) atanh(real(r)).*sqrt(d);
r2t  = @(r,d) real(r).*sqrt((d-2)./(1-real(r).^2));
%-----------------------------------

lw=1.2;
fs=15;

nRlz = 2000;
TVFlag = 'TVOff';
%TaperingMeth = 'ASAt_tc5';

T_list   = [100, 200, 600, 1200 1500 2000];
rho_list = [0,   0.2, 0.5, 0.7, 0.9];

%AC_list_lab={'W','AR1','AR4','AR14','AR-GPGR','AR-Yeo35-78','AR-Yeo85-100'};
AC_list_lab={'W','AR1','AR4','AR14','AR20A'};

TaperingMethList = {'v_ar1','v_ch','v_cr','v_ch4','v_cr4','ASAt','ASAt_ts1','ASAt_tt1','ASAt_tt2','ASAt_tc4','ASAt_tc5','v_roy'};

%TaperingMethList = {'ASAt_ts1'};


n_AC_list = numel(AC_list_lab);

for i = 1:n_AC_list
    for ii = 1:n_AC_list
        Labs{i,ii}=[AC_list_lab{i} ' - ' AC_list_lab{ii}];
    end
end

idx         = find(triu(ones(n_AC_list)));
idxdiag    = find(eye(n_AC_list));
idxoffdiag = setdiff(idx,idxdiag);

Mrk=[repmat('>ox',1,15)];
%repmat('o',1,11)]; 
%Mrk=Mrk(randperm(numel(idx_offdiag)));

spidx={[1 2],[3 4],[5 6],[7 8],[9 10],[11 12]};
for TaperingMeth=TaperingMethList
    
    fh_bias_SAS=figure('position',[50,500,810,350]);        hold on; %suptitle('Biases of Variance Estimators (ASAt)')
    % fh_asat_vars = figure('position',[50,500,1800,900]);      hold on; suptitle('Variances of ASAt')
    % 
    % fh_bias_CR=figure('position',[50,500,1800,900]);         hold on; suptitle('Biases of Variance Estimators (CR)')
    % 
    % fh_z_vars=figure('position',[50,500,1800,900]);          hold on; suptitle('Variance of Fisher Transformed Corr Coeffs')
    % fh_r_vars=figure('position',[50,500,1800,900]);          hold on; suptitle('Variance of Corr Coeffs')

    for t_cnt = 1:numel(T_list)    
        ndp = T_list(t_cnt);
        for rho_cnt = 1:numel(rho_list)
                rt_load='/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/R/';
                load([rt_load 'MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(rho_cnt).*100) '_' num2str(nRlz) '_xDF_' TVFlag '.mat'])
                
                %T_list   = [100, 200, 600, 1200];
                
                TarVar = eval(TaperingMeth{1});

                r0 = r0(1:n_AC_list,1:n_AC_list,:);
                TarVar = TarVar(1:n_AC_list,1:n_AC_list,:);
                
                %----VARIANCES-----------------------------------------------
                v_r0_tmp = var(r0,[],3);
                v_r0_D  (:,rho_cnt) = v_r0_tmp(idxdiag);
                v_r0_OD (:,rho_cnt) = v_r0_tmp(idxoffdiag);

                v_r0_Z_tmp  = var(atanh(r0),[],3);
                v_r0_Z_D (:,rho_cnt)  = v_r0_Z_tmp(idxdiag);
                v_r0_Z_OD (:,rho_cnt) = v_r0_Z_tmp(idxoffdiag);
                %----Biases of V----------------------------------------------
                v_ASAt_diag(:,rho_cnt)    = TarVar(idxdiag);
                v_ASAt_offdiag(:,rho_cnt) = TarVar(idxoffdiag);
                %--- SAS
                bv_SAS_tmp = real(bias(sqrt(var(r0,[],3)),sqrt(mean(TarVar,3))));
                bv_SAS_diag(:,rho_cnt)    = bv_SAS_tmp(idxdiag);
                bv_SAS_offdiag(:,rho_cnt) = bv_SAS_tmp(idxoffdiag);

                %--- CR
    %             bv_CR_tmp = bias(var(r0,[],3),mean(v_ch,3));
    %             bv_CR_diag(:,rho_cnt)    = bv_CR_tmp(idx_diag);
    %             bv_CR_offdiag(:,rho_cnt) = bv_CR_tmp(idx_offdiag);

                clear *_tmp
        end

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
    %     plot(((1-rho_list.^2).^2)./ndp,'linewidth',lw,'color',[.5 .5 .5],'marker','x') %theoritical variance
    %     
    %     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    %     ylabel ('$\sigma^2_r$','Interpreter','latex','fontsize',fs)
    %     
    % %     legend(Labs{idx_diag});    
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
    %     plot(((1-rho_list.^2).^2)./ndp,'linewidth',lw,'color',[.5 .5 .5],'marker','x') %theoritical variance
    %     
    %     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    %     ylabel ('$\sigma^2_r$','Interpreter','latex','fontsize',fs)

    %     legend(Labs{idx_offdiag});

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
    % %     legend(Labs{idx_diag});    
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

    %     legend(Labs{idx_offdiag});


        %------------------------------------------------------------------------------------
        %------------------------------------------------------------------------------------

    %     figure(fh_asat_vars);
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
    % %     legend(Labs{idx_diag});    
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

    %     legend(Labs{idx_offdiag});

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

        for i=1:numel(idxdiag)
            plot(bv_SAS_diag(i,:),'LineWidth',lw,'marker',Mrk(i));
        end

        if t_cnt>3
            xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
        end
        
        %xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
        %ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)

        %legend(Labs{idxdiag},'location','southwest');
        %------------------------------------------    
        sph11 = subplot(2,numel(T_list),spidx{t_cnt}(2)); 
        hold on; box on; grid on; 
        title(['T=' num2str(ndp) ', $\rho_{xx} \neq \rho_{yy}$' ],'Interpreter','latex')

        sph11.XTick      = 1:numel(rho_list);
        sph11.XTickLabel = {'0','.2','.5','.7','.9'};
        sph11.YTick = -50:20:50;
        ylim([-50 50]); xlim([1 numel(rho_list)]);
        for i=1:numel(idxoffdiag)
            plot(bv_SAS_offdiag(i,:),'LineWidth',lw,'marker',Mrk(i))
        end

        if t_cnt>3
            xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
        end
        
        %xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
        %ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)

        %legend(Labs{idxoffdiag},'location','southwest');
        %------------------------------------------------------------------------------------
        %------------------------------------------------------------------------------------

    %     figure(fh_bias_raw);
    %     sph10=subplot(2,4,spidx{t_cnt}(1)); 
    %     hold on; box on; grid on; 
    %     title(['T=' num2str(ndp) ', $\lambda_x=\lambda_y$' ],'Interpreter','latex')
    %     
    %     sph10.XTick      = 1:numel(rho_list);
    %     sph10.XTickLabel = {'0','.2','.5','.7','.9'};
    %     
    %     ylim([-150 40]); xlim([1 numel(rho_list)]);
    %     
    %     for i=1:numel(idx_diag)
    %         plot(bv_diag(i,:),'LineWidth',lw,'marker',Mrk(i));
    %     end
    %     
    %     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    %     ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    %     
    %     legend(Labs{idx_diag},'location','southwest');
    %     %------------------------------------------    
    %     sph11 = subplot(2,4,spidx{t_cnt}(2)); 
    %     hold on; box on; grid on; 
    %     title(['T=' num2str(ndp) ', $\lambda_x \neq \lambda_y$' ],'Interpreter','latex')
    %     
    %     sph11.XTick      = 1:numel(rho_list);
    %     sph11.XTickLabel = {'0','.2','.5','.7','.9'};
    %     
    %     ylim([-150 40]); xlim([1 numel(rho_list)]);
    %     for i=1:numel(idx_offdiag)
    %         plot(bv_offdiag(i,:),'LineWidth',lw,'marker',Mrk(i))
    %     end
    %     
    %     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    %     ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    %     
    %     legend(Labs{idx_offdiag},'location','southwest');
        %------------------------------------------------------------------------------------
        %------------------------------------------------------------------------------------

    %     figure(fh_bias_CR);
    %     sph10=subplot(2,4,spidx{t_cnt}(1)); 
    %     hold on; box on; grid on; 
    %     title(['T=' num2str(ndp) ', $\lambda_x=\lambda_y$' ],'Interpreter','latex')
    %     
    %     sph10.XTick      = 1:numel(rho_list);
    %     sph10.XTickLabel = {'0','.2','.5','.7','.9'};
    %     
    %     ylim([-100 100]); xlim([1 numel(rho_list)]);
    %     
    %     for i=1:numel(idx_diag)
    %         plot(bv_CR_diag(i,:),'LineWidth',lw,'marker',Mrk(i));
    %     end
    %     
    %     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    %     ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    %     
    % %     legend(Labs{idx_diag},'location','southwest');
    %     %------------------------------------------    
    %     sph11 = subplot(2,4,spidx{t_cnt}(2)); 
    %     hold on; box on; grid on; 
    %     title(['T=' num2str(ndp) ', $\lambda_x \neq \lambda_y$' ],'Interpreter','latex')
    %     
    %     sph11.XTick      = 1:numel(rho_list);
    %     sph11.XTickLabel = {'0','.2','.5','.7','.9'};
    %     
    %     ylim([-100 100]); xlim([1 numel(rho_list)]);
    %     for i=1:numel(idx_offdiag)
    %         plot(bv_CR_offdiag(i,:),'LineWidth',lw,'marker',Mrk(i))
    %     end
    %     
    %     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    %     ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    %     
    %     legend(Labs{idx_offdiag},'location','southwest');    

        %------------------------------------------------------------------------------------
        %------------------------------------------------------------------------------------

    %     figure(fh_bias_ar1);
    %     sph10=subplot(2,4,spidx{t_cnt}(1)); 
    %     hold on; box on; grid on; 
    %     title(['T=' num2str(ndp) ', $\lambda_x=\lambda_y$' ],'Interpreter','latex')
    %     
    %     sph10.XTick      = 1:numel(rho_list);
    %     sph10.XTickLabel = {'0','.2','.5','.7','.9'};
    %     
    %     ylim([-100 100]); xlim([1 numel(rho_list)]);
    %     
    %     for i=1:numel(idx_diag)
    %         plot(bv_ar1_diag(i,:),'LineWidth',lw,'marker',Mrk(i));
    %     end
    %     
    %     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    %     ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    %     
    %     legend(Labs{idx_diag},'location','southwest');
    %     %------------------------------------------    
    %     sph11 = subplot(2,4,spidx{t_cnt}(2)); 
    %     hold on; box on; grid on; 
    %     title(['T=' num2str(ndp) ', $\lambda_x \neq \lambda_y$' ],'Interpreter','latex')
    %     
    %     sph11.XTick      = 1:numel(rho_list);
    %     sph11.XTickLabel = {'0','.2','.5','.7','.9'};
    %     
    %     ylim([-100 100]); xlim([1 numel(rho_list)]);
    %     for i=1:numel(idx_offdiag)
    %         plot(bv_ar1_offdiag(i,:),'LineWidth',lw,'marker',Mrk(i))
    %     end
    %     
    %     xlabel('Correlation ($\rho_{xy}$)','Interpreter','latex')
    %     ylabel ('bias (\%)','Interpreter','latex','fontsize',fs)
    %     
    %     legend(Labs{idx_offdiag},'location','southwest');    
    end

    % set(fh_r_vars,'Color','w'); 
    % export_fig(fh_r_vars, 'Figs/MC/ME_r_vars.pdf')
    % 
    % set(fh_z_vars,'Color','w'); 
    % export_fig(fh_z_vars, 'Figs/MC/ME_z_vars.pdf')
    % 
    % set(fh_asat_vars,'Color','w'); 
    % export_fig(fh_asat_vars,'Figs/MC/ME_ASAt.pdf')
    % 
    %set(fh_bias_CR,'Color','w'); 
    %export_fig(fh_bias_CR,'Figs/MC/ME_bias_CnR.pdf')
    % 
    
    
    set(fh_bias_SAS,'Color','w'); 
    export_fig(fh_bias_SAS,['Figs/ME_bias_' TaperingMeth{1} '_' TVFlag '_' num2str(nRlz) '.pdf'])

    % set(fh_r_vars,'Color','w'); export_fig(fh_r_vars, 'Figs/MC/r_vars.pdf')
    % set(fh_z_vars,'Color','w'); export_fig(fh_z_vars, 'Figs/MC/z_vars.pdf')
    % 
    % set(fh_bias_ar1,'Color','w'); export_fig(fh_bias_ar1,'Figs/MC/bias_ar1.pdf')
    % set(fh_bias_raw,'Color','w'); export_fig(fh_bias_raw,'Figs/MC/bias_bcf_raw.pdf')
    % set(fh_bias_wd,'Color','w');  export_fig(fh_bias_wd, 'Figs/MC/bias_bcf_wd.pdf')
    % set(fh_bias_wdt,'Color','w'); export_fig(fh_bias_wdt,'Figs/MC/bias_bcf_wdt.pdf')
    
    %clear *_diag *_offdiag 
end
