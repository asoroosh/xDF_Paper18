clear

SiteList  = {'HCP'};
AtlasList = {'Gordon','Yeo','Power','ICA200'};
PP        = 'FPP';
save_rt   = '/Users/sorooshafyouni/Home/BCF/BCFAnal/Scrambling/HCP/ME/R/';

Col=get(groot,'defaultAxesColorOrder');

p2zfcn = @(p)-sqrt(2)*erfcinv(p*2);


alf_lv = [.005 0.025 0.05];
qz = -p2zfcn(alf_lv);

fs = 12;

nrndn   = 40;
subcomb = 325; 

idx = find(tril(ones(nrndn))); % diagnonal should be there?! what is the diagonal?

method_list       = {'naive','MEtt2','CR','ch' ,'AR1'};
method_list_label = {'Naive','xDF'  ,'BH','Q47','B35'};

n_mth = numel(method_list_label);

mthdCnt = 0.1;

for s=SiteList
    for a=AtlasList
        
        Histsfh = figure('position',[50,500,370,200]); hold on; box on; grid on; 
        QQsfh   = figure('position',[50,500,250,450]); hold on; box on; grid on; 
%         for cnt=1:subcomb
% 
%             rt = ['/Users/sorooshafyouni/Home/BCF/BCFAnal/Scrambling/HCP/ME/R/R_' s{1} '_' a{1} '_' PP '/'];
%             S  = load([rt 'InterSub_Scrambled_' s{1} '_' a{1} '_' PP '_' num2str(cnt) '.mat']);
% 
%             Z_naive = [Z_naive; S.Z_naive(idx)];
%             
%             Z_ME    = [Z_ME;    S.Z_ME(idx)];
%             Z_MEtt1 = [Z_MEtt1; S.Z_MEtt1(idx)];
%             Z_MEtt2 = [Z_MEtt2; S.Z_MEtt2(idx)];
%             Z_MEts1 = [Z_MEts1; S.Z_MEts1(idx)];
%             Z_MEtc5 = [Z_MEtc5; S.Z_MEtc5(idx)];
%             
%             Z_AR1   = [Z_AR1; S.Z_AR1(idx)];
%         end
%         
%         AllZs = [Z_naive Z_ME Z_MEtt1 Z_MEtt2 Z_MEts1 Z_MEtc5 Z_AR1];
        
        m_cnt = 1;
        for mth = method_list %because there are seven methods,
        %VarUnderHammer = AllZs(:,m_cnt);
        
            VarUnderHammer = [];
            for cnt=1:subcomb

                rt = ['/Users/sorooshafyouni/Home/BCF/BCFAnal/Scrambling/HCP/ME/R/R_' s{1} '_' a{1} '_' PP '/'];
                
                S  = load([rt 'InterSub_Scrambled_' s{1} '_' a{1} '_' PP '_' num2str(cnt) '_TVOff.mat']);

                VarNow = eval(['S.Z_' mth{1}]);
                VarUnderHammer = [VarUnderHammer; VarNow(idx)];
            end
        
            [~,~,KS(m_cnt)] = kstest(VarUnderHammer);
        
        %FPR---------------------------------------------------------------
            [mhm,XX] = HistPlotX_UniVarRef(VarUnderHammer,50,'figure',Histsfh,'color',Col(m_cnt,:));

    %         figure(fh); hold on; 
    %         
    %         line([q1z q1z] ,[0 mhm] ,'color',[.5 .5 .5 .5],'LineWidth',1.3)
    %         line(-[q1z q1z],[0 mhm] ,'color',[.5 .5 .5 .5],'LineWidth',1.3)
    %         line([q2z q2z] ,[0 mhm] ,'color',[.5 .5 .5 .5],'LineWidth',1.3)
    %         line(-[q2z q2z],[0 mhm] ,'color',[.5 .5 .5 .5],'LineWidth',1.3)
    %         
             for j = 1:3
                tq(j,m_cnt) = sum(abs(VarUnderHammer)>qz(j))./numel(VarUnderHammer);
             end
    %         scatter([tq1_tmp tq1_tmp] ,[mthdCnt mthdCnt] ,'markerfacecolor','b','markeredgecolor','b','marker','x','LineWidth',1.3);
    %         scatter(-[tq1_tmp tq1_tmp],[mthdCnt mthdCnt],'markerfacecolor','b','markeredgecolor','b','marker','x','LineWidth',1.3);
    %         scatter([tq2_tmp tq2_tmp] ,[mthdCnt mthdCnt] ,'markerfacecolor','r','markeredgecolor','r','marker','x','LineWidth',1.3);
    %         scatter(-[tq2_tmp tq2_tmp],[mthdCnt mthdCnt],'markerfacecolor','r','markeredgecolor','r','marker','x','LineWidth',1.3);
        %------------------------------------------------------------------

        %QQplots-----------------------------------------------------------
            figure(QQsfh)
            hqq=qqplot(VarUnderHammer([1:50 51:50:length(VarUnderHammer)-50 length(VarUnderHammer)-49:length(VarUnderHammer)]));
            hqq(1).Marker='.'; 
            hqq(1).MarkerSize=8;
            hqq(1).MarkerEdgeColor=Col(m_cnt,:);
            hqq(3).Color=Col(m_cnt,:); 
            hqq(3).LineStyle='-.';
            
            m_cnt = m_cnt+1;
        end        
        
        figure(Histsfh)
        %legend(method_list_label,'fontsize',fs)
        %plot(XX,normpdf(XX),'color',[1 0 0 0.4],'linewidth',2.7,'linestyle','-.')
        xlim([-15 15])
        xlabel('z-statistics','fontsize',fs)
        ylabel('Probability','fontsize',fs)
        legend(method_list_label,'fontsize',fs)
        
        figure(QQsfh)
        rlh = refline(1,0);
        rlh.Color = 'k'; 
        rlh.LineWidth = 1.3;
        rlh.LineStyle = '-.';
        title(['']);
        %xlim([-7 7])
        %ylim([-7 7])
        %legend(method_list_label,'fontsize',fs)
        
        disp([s{1} '_' a{1}])
        %save([save_rt 'BtwnSubBtwnNd_ME_' s{1} '_' a{1} '.mat'],'Z_naive','Z_MEts1');
        clear Z_naive Z_MEts1 *_tmp
        
        
        
        set(Histsfh,'Color','w');

        % 

        figure(QQsfh)
        ylabel('z-statistics','fontsize',fs,'interpreter','latex')
        xlabel('Standard Normal','fontsize',fs,'interpreter','latex')
        set(QQsfh,'Color','w');

        %---KS Statistics
        bfh = figure('position',[50,500,200,200]); 
        hold on; box on; 
        for i=1:numel(method_list_label)
            KSlog = -log10(KS(i));
            barfh = bar(i,KSlog,'FaceColor',Col(i,:),'FaceAlpha',0.5);
            txth = text(i,KSlog+0.1,num2str(round(KSlog,2)));
            set(txth,'rotation',90)
        end
        bfh.Children.XTick=[1:numel(method_list_label)];
        bfh.Children.XTickLabel=method_list_label; %MCPS: Monte Carlo Parametric Simulations
        bfh.Children.XTickLabelRotation=45;
        %bfh.Children.XTick = [1:3];
        ylabel('-log10(KS Statistics)','fontsize',fs,'interpreter','latex')
        xlim([0 n_mth+1])
        ylim([0 3.5])
        %barfh.FaceAlpha=0.3;
        %barfh.FaceColor='k';
        %barfh.EdgeColor='k';
        set(bfh,'Color','w');

        %---FPR Analysis

        fh_fpr = figure('position',[50,500,520,270]); 
        for j = 1:3
            sph_fpr = subplot(1,3,j);
            hold on; box on;
            for i=1:numel(method_list_label)
                tq_tmp = tq(j,i);
                barfh = bar(i,tq(j,i),'FaceColor',Col(i,:),'FaceAlpha',0.5);
                txth = text(i,tq_tmp+0.01,num2str(round(tq_tmp,3)));
                set(txth,'rotation',90)
            end
            line([0 n_mth+1],[alf_lv(j) alf_lv(j)].*2,'color','r','linestyle','-.')
            sph_fpr.XTick              = [1:numel(method_list_label)];
            sph_fpr.XTickLabel         = method_list_label; %MCPS: Monte Carlo Parametric Simulations
            sph_fpr.XTickLabelRotation = 45;
            %bfh.Children.XTick = [1:3];
            ylabel('False Positive Rate','fontsize',fs,'interpreter','latex')
            xlim([0 n_mth+1])
            ylim([0 0.4])
            %barfh.FaceAlpha=0.3;
            %barfh.FaceColor='k';
            %barfh.EdgeColor='k';
        end
        set(fh_fpr,'Color','w');


        fh_fpr2 = figure('position',[50,500,200,200]); 
        hold on; box on;
        for i=1:numel(method_list_label)
            tq_tmp = tq(2,i);
            barfh = bar(i,tq_tmp,'FaceColor',Col(i,:),'FaceAlpha',0.5);
            txth = text(i,tq_tmp+0.01,num2str(round(tq_tmp,3)));
            set(txth,'rotation',90)
        end
        line([0 n_mth+1],[alf_lv(2) alf_lv(2)].*2,'color','r','linestyle','-.')
        fh_fpr2.Children.XTick              = [1:numel(method_list_label)];
        fh_fpr2.Children.XTickLabel         = method_list_label; %MCPS: Monte Carlo Parametric Simulations
        fh_fpr2.Children.XTickLabelRotation = 45;
        %bfh.Children.XTick = [1:3];
        ylabel('False Positive Rate','fontsize',fs,'interpreter','latex')
        xlim([0 n_mth+1])
        ylim([0 0.4])
        %barfh.FaceAlpha=0.3;
        %barfh.FaceColor='k';
        %barfh.EdgeColor='k';
        set(fh_fpr2,'Color','w');

        %%%%EXPORT
        export_fig(QQsfh,['Figs/' s{1} '_' a{1} '_QQ.pdf'])
        export_fig(bfh,['Figs/' s{1} '_' a{1} '_KS.pdf'])
        export_fig(Histsfh,['Figs/' s{1} '_' a{1} '_Hist.pdf'])
        export_fig(fh_fpr2,['Figs/' s{1} '_' a{1} '_FPR5p.pdf'])
        export_fig(fh_fpr,['Figs/' s{1} '_' a{1} '_FPR.pdf'])
    end
end

% 

% bfh_fpr = figure('position',[50,500,200,200]); 
% hold on; box on;
% for i=1:numel(method_list_label)
%     barfh = bar(i,tq2(i),'FaceColor',Col(i,:),'FaceAlpha',0.5);
% end
% line([0 n_mth+1],[.05 .05],'color','r','linestyle','-.')
% bfh_fpr.Children.XTick=[1:numel(method_list_label)];
% bfh_fpr.Children.XTickLabel=method_list_label; %MCPS: Monte Carlo Parametric Simulations
% bfh_fpr.Children.XTickLabelRotation=45;
% %bfh.Children.XTick = [1:3];
% ylabel('False Positive Rate','fontsize',fs,'interpreter','latex')
% xlim([0 n_mth+1])
% %barfh.FaceAlpha=0.3;
% %barfh.FaceColor='k';
% %barfh.EdgeColor='k';
% set(bfh_fpr,'Color','w');
% 
% bfh_fpr = figure('position',[50,500,200,200]); 
% hold on; box on;
% for i=1:numel(method_list_label)
%     barfh = bar(i,tq3(i),'FaceColor',Col(i,:),'FaceAlpha',0.5);
% end
% line([0 n_mth+1],[.10 .10],'color','r','linestyle','-.')
% bfh_fpr.Children.XTick=[1:numel(method_list_label)];
% bfh_fpr.Children.XTickLabel=method_list_label; %MCPS: Monte Carlo Parametric Simulations
% bfh_fpr.Children.XTickLabelRotation=45;
% %bfh.Children.XTick = [1:3];
% ylabel('False Positive Rate','fontsize',fs,'interpreter','latex')
% xlim([0 n_mth+1])
% %barfh.FaceAlpha=0.3;
% %barfh.FaceColor='k';
% %barfh.EdgeColor='k';
% set(bfh_fpr,'Color','w');
%-----

% export_fig(bfh_fpr,['FPR.pdf'])
% 
% export_fig(bfh,['KS.pdf'])
% 
% set(Histsfh,'Color','w');
% export_fig(Histsfh,['Hist.pdf'])
% 
% set(QQsfh,'Color','w');
% export_fig(QQsfh,['QQ.pdf'])    
% end
                

