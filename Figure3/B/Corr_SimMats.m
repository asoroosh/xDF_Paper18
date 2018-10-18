clear

nn=114; 

lw = 1.3;
fs = 12; 

TVFlag = 'TVOff';

%SubID = 118528; 
SubID = 135932; 

addpath('/Users/sorooshafyouni/Home/matlab/Ext/contest')

Col = get(groot,'defaultAxesColorOrder');


%% Graph
DensRng=0.01:0.01:0.50;

nRlz = 2000;

idx = find(triu(ones(nn),1));

for rlz=1:nRlz
    if ~mod(rlz,100); disp(num2str(rlz)); end;
    load(['R/RndNetworks_Yeo_FPP_' num2str(SubID) '_' num2str(rlz) '_' TVFlag '_Mats.mat'])
        %Matrices
        All_adjmats_Naive (:,rlz)       = real(Mat.amat(idx));
        All_rmats (:,rlz)       = real(Mat.rmat(idx));
        All_adjmats_ME (:,rlz)  = real(Mat.adjt2_amat(idx)); 
        All_adjmats_CR (:,rlz)  = real(Mat.adjCR_amat(idx));
        All_adjmats_CH (:,rlz)  = real(Mat.adjCH_amat(idx));
        All_adjmats_Fox (:,rlz) = real(Mat.adjFox_amat(idx));
        All_adjmats_FSL (:,rlz) = real(Mat.adjFSL_amat(idx));
end
m_All_adjamats_FSL = All_adjmats_FSL(:);
m_All_adjamats_Fox = All_adjmats_Fox(:);
m_All_adjamats_CR = All_adjmats_CR(:);
m_All_adjamats_CH = All_adjmats_CH(:);
m_All_adjamats_ME = All_adjmats_ME(:);
m_All_adjamats_Naive = All_adjmats_Naive(:);
m_All_rmats = All_rmats(:);
%--------------------------------------------------------------------------

method_label    = {'Naive','xDF','BH','Q47','AR1MCPS','G-Q47'};
method_list     = {'Naive','ME','CR','CH','FSL','Fox'};

n_mth = numel(method_label);

p2zfcn = @(p)-sqrt(2)*erfcinv(p*2);

alf_lv = [0.005 0.025 0.05];
qz = -p2zfcn(alf_lv);

for i = 1:3
    m_cnt = 1;
    for mth = method_list
        MthVar = eval(['m_All_adjamats_' mth{1}]);
        tq(i,m_cnt) = sum(abs(MthVar) >qz(i))./numel(MthVar); 
        m_cnt = m_cnt + 1;
    end
end


fh_fpr = figure('position',[50,500,520,270]); 
for j = 1:3
    sph_fpr = subplot(1,3,j);
    hold on; box on;
    for i=1:n_mth
        tq_tmp = tq(j,i);
        barfh = bar(i,tq_tmp,'FaceColor',Col(i,:),'FaceAlpha',0.5);
        
        txth = text(i,tq_tmp+0.01,num2str(round(tq_tmp,3)));
        set(txth,'rotation',90)
    end
    line([0 n_mth+1],[alf_lv(j) alf_lv(j)].*2,'color','r','linestyle','-.')
    sph_fpr.XTick              = [1:n_mth];
    sph_fpr.XTickLabel         = method_label; %MCPS: Monte Carlo Parametric Simulations
    sph_fpr.XTickLabelRotation = 45;
    %bfh.Children.XTick = [1:3];
    ylabel('False Positive Rate','fontsize',fs,'interpreter','latex')
    xlim([0 n_mth+1])
    ylim([0 0.50])
end
set(fh_fpr,'Color','w');

fh_fpr2 = figure('position',[50,500,200,200]); 
hold on; box on;
for i=1:numel(method_label)
    tq_tmp = tq(2,i);
    barfh = bar(i,tq_tmp,'FaceColor',Col(i,:),'FaceAlpha',0.5);
    txth = text(i,tq_tmp+0.01,num2str(round(tq_tmp,3)));
    set(txth,'rotation',90)
    
end
line([0 n_mth+1],[alf_lv(2) alf_lv(2)].*2,'color','r','linestyle','-.')
fh_fpr2.Children.XTick              = [1:numel(method_label)];
fh_fpr2.Children.XTickLabel         = method_label; %MCPS: Monte Carlo Parametric Simulations
fh_fpr2.Children.XTickLabelRotation = 45;
%bfh.Children.XTick = [1:3];
ylabel('False Positive Rate','fontsize',fs,'interpreter','latex')
xlim([0 n_mth+1])
%barfh.FaceAlpha=0.3;
%barfh.FaceColor='k';
%barfh.EdgeColor='k';
ylim([0 0.50])
set(fh_fpr2,'Color','w');


%--------------------------------------------------------------------------
fh_qq = figure('position',[50,500,250,450]);
hold on; grid on; box on;

smplstp = 1000;

hqq_ac=qqplot(m_All_adjamats_Naive(1:smplstp:end));
hqq_ac(1).Marker='.'; hqq_ac(1).MarkerEdgeColor=Col(1,:); hqq_ac(3).Color='none';

hqq_adj_me=qqplot(m_All_adjamats_ME(1:smplstp:end));
hqq_adj_me(1).Marker='.'; hqq_adj_me(1).MarkerEdgeColor=Col(2,:); hqq_adj_me(3).Color='none';

hqq_adj_cr=qqplot(m_All_adjamats_CR(1:smplstp:end));
hqq_adj_cr(1).Marker='.'; hqq_adj_cr(1).MarkerEdgeColor=Col(3,:); hqq_adj_cr(3).Color='none';

hqq_adj_ch=qqplot(m_All_adjamats_CH(1:smplstp:end));
hqq_adj_ch(1).Marker='.'; hqq_adj_ch(1).MarkerEdgeColor=Col(4,:); hqq_adj_ch(3).Color='none';

hqq_adj_fox=qqplot(m_All_adjamats_Fox(1:smplstp:end));
hqq_adj_fox(1).Marker='.'; hqq_adj_fox(1).MarkerEdgeColor=Col(6,:); hqq_adj_fox(3).Color='none';

hqq_adj_fsl=qqplot(m_All_adjamats_FSL(1:smplstp:end));
hqq_adj_fsl(1).Marker='.'; hqq_adj_fsl(1).MarkerEdgeColor=Col(5,:); hqq_adj_fsl(3).Color='none';

title('')

xlim([-5 5])
ylim([-5 5])

%legend([hqq_ac(1), hqq_adj_cr(1), hqq_adj_ch(1), hqq_adj_me(1), hqq_adj_fox(1), hqq_adj_fsl(1)],...
    %{'Naive','BH','B46','xDF','G-B46','AR1MCPS'},'fontsize',fs,'location','southeast')

plot(-7:0.5:7,-7:0.5:7,'k-.')

ylabel('z-statistics','fontsize',fs,'interpreter','latex')
xlabel('Standard Normal','fontsize',fs,'interpreter','latex')
set(gcf,'Color','w');
%--------------------------------------------------------------------------
bfh_hist = figure('position',[50,500,370,200]); 
hold on; grid on; box on;

m_cnt = 1;
for mth = method_list
    MthVar = eval(['m_All_adjamats_' mth{1}]);
    HistPlotX_UniVarRef(MthVar,100,'figure',bfh_hist,'color',Col(m_cnt,:));
    m_cnt = m_cnt + 1; 
end
legend(method_label,'fontsize',fs)
set(gcf,'Color','w');
%--------------------------------------------------------------------------
[~,~,KS_FSL]   = kstest(All_adjmats_FSL);
[~,~,KS_Fox]   = kstest(All_adjmats_Fox);
[~,~,KS_CR]    = kstest(All_adjmats_CR);
[~,~,KS_CH]    = kstest(All_adjmats_CH); 
[~,~,KS_ME]    = kstest(All_adjmats_ME); 
[~,~,KS_Naive] = kstest(All_adjmats_Naive);
%[~,~,KS_rnd] = kstest(All_rmats);
% 
bfh_KS = figure('position',[50,500,200,200]); 
%bfh_KS =figure('position',[50,500,220,270]);
hold on; box on;

m_cnt = 1; 
for mth = method_list
    MthVar = eval(['All_adjmats_' mth{1}]);
    [~,~,KS] = kstest(MthVar);
    KSlog = -log10(KS);
    barfh = bar(m_cnt,KSlog,'FaceColor',Col(m_cnt,:),'FaceAlpha',0.5);
    
    txth = text(m_cnt,KSlog+0.1,num2str(round(KSlog,2)));
    set(txth,'rotation',90)
    
    m_cnt = m_cnt+1;
end
bfh_KS.Children.XTick=[1:numel(method_label)];
bfh_KS.Children.XTickLabel=method_label; %MCPS: Monte Carlo Parametric Simulations
bfh_KS.Children.XTickLabelRotation=45;
ylabel('-$\log$10(KS Statistics)','fontsize',fs,'interpreter','latex')
xlim([0 numel(method_list)+1])
set(bfh_KS,'Color','w');
ylim([0 3.5])

%-----------EXPORTS
% export_fig(fh_gcrate,'Topo.pdf')

export_fig(fh_fpr,  ['FPR_' num2str(SubID) '.pdf']);
export_fig(fh_fpr2, ['FPR5p_' num2str(SubID) '.pdf']);
export_fig(bfh_hist,['Hist_' num2str(SubID) '.pdf']);
export_fig(bfh_KS,  ['KS_' num2str(SubID) '.pdf']);
export_fig(fh_qq,   ['QQ_' num2str(SubID) '.pdf']);

% export_fig(fh_gcrate,'GChanges.pdf');