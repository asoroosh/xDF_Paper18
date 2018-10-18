clear
%close all

% We use to add two options here but we later removed them. 
%
% 1) Non parametric testings, following the model in the paper; I removed
% it because it is not the focus. The results are great, but it is because
% we follow the model which has AC structures in it. 
%
% 2) Parametric Monte-carlo, I also removed it, because it is not
% computationally efficient. Note, this option is what FSLnets does, BUT
% with much crappier model (AR1). 
%   
%---------ANNON FUNCs---------------
% i.e. negative bias means UNDER estimated
%     positive bias means OVER estimated
% E : estimated (somehow)
% C : Calculated (i.e. your ground-truth)
bias = @(MC_C,ParamE) (ParamE-MC_C)./MC_C.*100;

%bias = @(old,new) (new-old)./old.*100; %Tom's version!
r2z  = @(r,d) atanh(real(r)).*sqrt(d);
r2t  = @(r,d) real(r).*sqrt((d-2)./(1-real(r).^2));
%-----------------------------------

% alp = 0.05;
% alp = alp/2;

alpList = [0.01 0.05 0.10];

lw = 1.2;
fs = 12;

nRlz = 2000;
TVFlag = 'TVOff';

T_list   = [100, 200, 600, 1200, 1500, 2000];
rho_list = [0, 0.2, 0.5, 0.7, 0.9];

AC_list_lab = {'W','AR1','AR4','AR14','AR20A'};

TaperingMethList =  {'ar1','ch' ,'cr','f'  ,'f_ts1','f_tt1','f_tt2','f_tc4','f_tc5'};
%TaperingMethList =  {'f_tt1'};
MethodNames =       {'B35','Q47','BH','xDF','xDF'  ,'xDF'  ,'xDF','xDF','xDF'};

%TaperingMethList = {'f'};
%MethodNames = {'xDFNT'};

n_AC_list = numel(AC_list_lab);

for i = 1:n_AC_list
    for ii = 1:n_AC_list
        Labs{i,ii}=[AC_list_lab{i} ' - ' AC_list_lab{ii}];
    end
end

idx         = find(triu(ones(n_AC_list)));
idx_diag    = find(eye(n_AC_list));
idx_offdiag = setdiff(idx,idx_diag);

%Mrk=[repmat('>o',1,11)]; %Mrk=Mrk(randperm(numel(idx_offdiag)));
Mrk=[repmat('>ox',1,15)]; 

spidx={[1 2],[3 4],[5 6],[7 8]};

alp_cnt = 1; 
for alp = alpList
    tm_cnt = 1 ;
    for TaperingMeth=TaperingMethList

        for t_cnt = 1:numel(T_list)    
            %for rho_cnt = 1:numel(rho_list)
                rt_load = '/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/R/';
                load([rt_load 'MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(1).*100) '_' num2str(nRlz) '_xDF_' TVFlag '.mat'])

                %T_list   = [100, 200, 600, 1200, 1500 2000];
                ndp = T_list(t_cnt);
                
                if strcmp(TaperingMeth{1},'f')
                    ['p_' TaperingMeth{1}]
                    TarVar = eval(['p_' TaperingMeth{1}]);
                    TarVar = TarVar(1:n_AC_list,1:n_AC_list,:);
                    fpr_p_f_tmp = sum(TarVar<alp,3)./nRlz;
                else
                    TarVar = eval(['z_' TaperingMeth{1}]);
                    TarVar = TarVar(1:n_AC_list,1:n_AC_list,:);
                    fpr_p_f_tmp = sum(2*normcdf(-abs(TarVar))<alp,3)./nRlz;
                end
                
                r0 = r0(1:n_AC_list,1:n_AC_list,:);
                
                 %fpr_p_f_tmp = sum(TarVar<alp,3)./nRlz;
                fpr_p_f_d (:,t_cnt) = fpr_p_f_tmp(idx_diag)*100;
                fpr_p_f_od(:,t_cnt) = fpr_p_f_tmp(idx_offdiag)*100;
                
                zn = atanh(r0).*sqrt(ndp-3);
                fpr_pvn_tmp = sum(2*normcdf(-abs(zn))<alp,3)./nRlz;
                fpr_pvn_d (:,t_cnt) = fpr_pvn_tmp(idx_diag)*100;
                fpr_pvn_od(:,t_cnt) = fpr_pvn_tmp(idx_offdiag)*100;

        %         load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FPR/R/FPR_MESim_t' num2str(T_list(t_cnt)) '_r' num2str(rho_list(1).*100) '.mat'],'p_np_f_ts1','p_np_fisher')
        %         fpr_np_f_tmp = sum(p_np_f_ts1<0.05,3)./nRlz;        
        %         fpr_np_f_d (:,t_cnt) = fpr_np_f_tmp(idx_diag)*100;
        %         fpr_np_f_od(:,t_cnt)= fpr_np_f_tmp(idx_offdiag)*100;
        %         
        %         fpr_np_fisher_tmp = sum(p_np_fisher<0.05,3)./nRlz;        
        %         fpr_np_fisher_d (:,t_cnt) = fpr_np_fisher_tmp(idx_diag)*100;
        %         fpr_np_fisher_od(:,t_cnt) = fpr_np_fisher_tmp(idx_offdiag)*100;
            %end 
        end

        FPR_D{alp_cnt,tm_cnt}  = fpr_p_f_d;
        FPR_OD{alp_cnt,tm_cnt} = fpr_p_f_od;
        
        %-------------------------------------------
        %-------------------------------------------
        FPRfh = figure('position',[50,500,380,350]);        
        hold on; %suptitle('Specificity');

        figure(FPRfh);
        sph00=subplot(2,2,1);
        %subplot(2,4,spidx{t_cnt}(1)); 
        hold on; box on; grid on; 
        title(['Naive - $\rho_{xx} = \rho_{yy}$'],'Interpreter','latex')

        sph00.XTick      = 1:numel(T_list);
        sph00.XTickLabel = {'100','200','600','1200','1500','2000'};

        sph00.YTick      = [1 5 10 20:20:100];
        sph00.YTickLabel = num2str([1 5 10 20:20:100]');

        ylim([0 60]); 
        xlim([1 numel(T_list)]);

        for i=1:numel(idx_diag)
            plot(fpr_pvn_d(i,:),'LineWidth',lw,'marker',Mrk(i));
        end

        xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
        ylabel ('FP Rate \%','Interpreter','latex','fontsize',fs)

        %legend(Labs{idx_diag});    
        %-----------------------------------
        % sph00=subplot(2,4,2);
        % %subplot(2,4,spidx{t_cnt}(1)); 
        % hold on; box on; grid on; 
        % title(['Fisher Non-Parametric- $\lambda_x = \lambda_y$'],'Interpreter','latex')
        % 
        % sph00.XTick      = 1:numel(T_list);
        % sph00.XTickLabel = {'100','200','600','1200'};
        % 
        % ylim([0 50]); 
        % xlim([1 numel(T_list)]);
        % 
        % for i=1:numel(idx_diag)
        %     plot(fpr_np_fisher_d(i,:),'LineWidth',lw,'marker',Mrk(i));
        % end
        % 
        % xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
        % ylabel ('FP Rate \%','Interpreter','latex','fontsize',fs)

        %legend(Labs{idx_diag});    
        %-----------------------------------
        sph02 = subplot(2,2,2);
        hold on; box on; grid on; 
        title([MethodNames{tm_cnt} ' - $\rho_{xx} = \rho_{yy}$' ],'Interpreter','latex')

        sph02.YTick      = [1 5 10 20:20:100];
        sph02.YTickLabel = num2str([1 5 10 20:20:100]');

        sph02.XTick      = 1:numel(T_list);
        sph02.XTickLabel = {'100','200','600','1200','1500','2000'};

        ylim([0 20]); 
        xlim([1 numel(T_list)]);

        for i=1:numel(idx_diag)
            plot(fpr_p_f_d(i,:),'LineWidth',lw,'marker',Mrk(i));
        end

        xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
        ylabel ('FP Rate \%','Interpreter','latex','fontsize',fs)

        %legend(Labs{idx_diag});    
        %-----------------------------------
        % sph02 = subplot(2,4,4);
        % hold on; box on; grid on; 
        % title(['Monster Non-Parametric - $\lambda_x = \lambda_y$' ],'Interpreter','latex')
        % 
        % sph02.XTick      = 1:numel(T_list);
        % sph02.XTickLabel = {'100','200','600','1200'};
        % 
        % ylim([0 50]); 
        % xlim([1 numel(T_list)]);
        % 
        % for i=1:numel(idx_diag)
        %     plot(fpr_np_f_d(i,:),'LineWidth',lw,'marker',Mrk(i));
        % end
        % 
        % xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
        % ylabel ('FP Rate \%','Interpreter','latex','fontsize',fs)
        % 
        % legend(Labs{idx_diag});    

        %--------------------------------------------------------------------------  
        %--------------------------------------------------------------------------

        sph01 = subplot(2,2,3);
        hold on; box on; grid on; 
        title(['Naive - $\rho_{xx} \neq \rho_{yy}$'],'Interpreter','latex')

        sph01.XTick      = 1:numel(T_list);
        sph01.XTickLabel = {'100','200','600','1200','1500','2000'};

        sph01.YTick      = [1 5 10 20:20:100];
        sph01.YTickLabel = num2str([1 5 10 20:20:100]');

        ylim([0 60]); 
        xlim([1 numel(T_list)]);

        for i=1:numel(idx_offdiag)
            plot(fpr_pvn_od(i,:),'LineWidth',lw,'marker',Mrk(i));
        end

        xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
        ylabel ('FP Rate \%','Interpreter','latex','fontsize',fs)

        %legend(Labs{idx_offdiag});

        %------------------------------------------    

        % sph01 = subplot(2,4,6);
        % hold on; box on; grid on; 
        % title(['Fisher Non-Parametric - $\lambda_x \neq \lambda_y$'],'Interpreter','latex')
        % 
        % sph01.XTick      = 1:numel(T_list);
        % sph01.XTickLabel = {'100','200','600','1200'};
        % 
        % ylim([0 50]); 
        % xlim([1 numel(T_list)]);
        % 
        % for i=1:numel(idx_offdiag)
        %     plot(fpr_np_fisher_od(i,:),'LineWidth',lw,'marker',Mrk(i));
        % end
        % 
        % xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
        % ylabel ('FP Rate \%','Interpreter','latex','fontsize',fs)

        %legend(Labs{idx_offdiag});

        %--------------------------------------------

        sph03 = subplot(2,2,4);
        hold on; box on; grid on; 
        title([MethodNames{tm_cnt} ' - $\rho_{xx} \neq \rho_{yy}$' ],'Interpreter','latex')

        sph03.YTick      = [1 5 10 20:20:100];
        sph03.YTickLabel = num2str([1 5 10 20:20:100]');

        sph03.XTick      = 1:numel(T_list);
        sph03.XTickLabel = {'100','200','600','1200','1500','2000'};

        ylim([0 20]); 
        xlim([1 numel(T_list)]);

        for i=1:numel(idx_offdiag)
            plot(fpr_p_f_od(i,:),'LineWidth',lw,'marker',Mrk(i));
        end

        xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
        ylabel ('FP Rate \%','Interpreter','latex','fontsize',fs)

        %legend(Labs{idx_offdiag});
        %------------------------------------------    

        % sph03 = subplot(2,4,8);
        % hold on; box on; grid on; 
        % title(['Monster Non-Parametric - $\lambda_x \neq \lambda_y$' ],'Interpreter','latex')
        % 
        % sph03.XTick      = 1:numel(T_list);
        % sph03.XTickLabel = {'100','200','600','1200'};
        % 
        % ylim([0 50]); 
        % xlim([1 numel(T_list)]);
        % 
        % for i=1:numel(idx_offdiag)
        %     plot(fpr_np_f_od(i,:),'LineWidth',lw,'marker',Mrk(i));
        % end
        % 
        % xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
        % ylabel ('FP Rate \%','Interpreter','latex','fontsize',fs)
        % 
        % legend(Labs{idx_offdiag});

        %------------------------------------------

        set(FPRfh,'Color','w');
        export_fig(FPRfh,['Figs/FPR_' TaperingMeth{1} '_' num2str(alp) '_' TVFlag '_' num2str(nRlz) '.pdf'])
        
        tm_cnt = tm_cnt+1;
        
    end
    alp_cnt = alp_cnt+1;
end

fh_mfprb_diag=figure('position',[50,500,600,250]);
for alp_cnt = 1:numel(alpList)
    subplot(1,3,alp_cnt); hold on; box on; grid on; 
    title(['$\alpha$-level = ' num2str(alpList(alp_cnt)*100) '\%'],'Interpreter','latex','FontSize',fs)
    for mth_cnt = [1 2 3]
        MeanFPRBias = mean(bias(ones(n_AC_list,numel(T_list)).*alpList(alp_cnt).*100,FPR_D{alp_cnt,mth_cnt}));
        plot(T_list,MeanFPRBias,'linewidth',lw,'marker','o')
       clear MeanFPRBias
    end
    
    line([0 1300],[0 0],'color','r','linestyle',':')
    
    ylabel('Mean Bias \% of FPR','Interpreter','latex','FontSize',fs)
    xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
    legend(MethodNames{[1 2 3]})
    ylim([-50 150])
    xlim([0 1300])
end
set(gcf,'color','w')
%export_fig(fh_mfprb_diag,'MeanFPRBias_Diag.pdf')

fh_mfprb_offdiag = figure('position',[50,500,600,250]);
for alp_cnt = 1:numel(alpList)
    subplot(1,3,alp_cnt); hold on; box on; grid on; 
    title(['$\alpha$-level = ' num2str(alpList(alp_cnt)*100) '\%'],'Interpreter','latex','FontSize',fs)    
    for mth_cnt = [1 2 3]
        MeanFPRBias = mean(bias(ones(10,numel(T_list)).*alpList(alp_cnt).*100,FPR_OD{alp_cnt,mth_cnt}));
        plot(T_list,MeanFPRBias,'linewidth',lw,'marker','o')
       clear MeanFPRBias
    end
    
    line([0 1300],[0 0],'color','r','linestyle',':')
    
    ylabel('Mean Bias \% of FPR','Interpreter','latex','FontSize',fs)
    xlabel('Length ($T$)','Interpreter','latex','FontSize',fs)
    legend(MethodNames{[1 2 3]})
    ylim([-50 150])
    xlim([0 1300])
end
set(gcf,'color','w')
%export_fig(fh_mfprb_offdiag,'MeanFPRBias_OffDiag.pdf')
