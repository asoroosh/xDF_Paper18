%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reproduces Figure 1: Inter-subject correlation between two HCP subjects
%
%%% REQUIREMENTS:
% 1) Data: HCP ICA-fixed volumes subject 135932 & 118528, parcellated with Yeo atlas
% 2) Code: xDF package, available via: https://github.com/asoroosh/xDF/
% 
% Soroosh Afyouni, University of Oxford, 2019, 
% srafyouni@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/xDF'))

fs=12;

lw = 1.3;

xxlim=7;

rightmxxlim = 50; 

z = @(p) -sqrt(2) * erfcinv(p*2);
alp=0.05;

p2z = -z(alp/2);

T = 1200;

TR = 0.72;
scanses = round(T*TR);
timspan = TR:TR:scanses;

Col=get(groot,'defaultAxesColorOrder');
Blue  = Col(1,:); % Blue
Yellow  = Col(3,:); % Yellow
Red     = Col(2,:); % Red (/orange!)
Purple  = Col(4,:); % Purple
Green   = Col(5,:); % Green
CRCol   = Col(6,:);

% load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/Pocket/RealData_ts_ac.mat')
% AdjCUR = UNRATE;

load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/HetBivStressTest/FC135932/HCP_FPP_118528_OnlyMTS.mat')
AdjCUR = mts(:,104); %104 in 118528 has the largest autocorrelation length! 
AdjCUR = (AdjCUR-mean(AdjCUR))./std(AdjCUR);
clear mts

load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/HetBivStressTest/FC135932/HCP_FPP_135932_Yeo_ROIs.mat','mts')

ph0 = figure('position',[50,500,1100,200]); 
hold on; box on; grid on; 
title('HCP118528 Left Dorsal Prefrontal Cortex (PFCd)','fontsize',fs);
plot(timspan,AdjCUR,'linewidth',lw,'marker','.','color','k');
%plot(CUR,'linewidth',lw,'marker','none','color',[.5 .5 .5 .5]);
ylabel('Standardised BOLD Response','fontsize',fs,'interpreter','latex')
xlabel('Seconds','fontsize',fs,'interpreter','latex')
%ph0.Children.XTick=1:24:840;
%ph0.Children.XTickLabel=num2cell(1949:2:2018);
%ph0.Children.XTickLabelRotation=45;
xlim([0 max(timspan)])
%legend({'UR(Seasonally Adjusted)','UR'},'fontsize',fs,'location','northwest')
set(ph0,'Color','w');

T=numel(AdjCUR);

[AC_RWV,bnd] = AC_fft(AdjCUR,T);
bnd = bnd(1);

ts=mts(1:T,:);

% If you try pre-whitening, then ALL points lay on the reference line...
%
% disp('Prewhitening the resting-state data...')
% ts = PreWhitenMe(ts,T,'taper','tukey',sqrt(T))';

[r,p]=corr(AdjCUR,ts);

sc_naive_rate = sum(fdr_bh(p))./114*100;

znaive = atanh(r).*sqrt(numel(AdjCUR)-3);
% histogram(znaive,114,'Normalization','probability')
% line([1.64 1.64],[0 .1])
% line(-[1.64 1.64],[0 .1])

ac_ts=AC_fft(ts,numel(AdjCUR));
ac_strength = sum((ac_ts(:,1:round(T/9)).*AC_RWV(1:round(T/9))).^2,2);

for i=1:114
    %[~,Stat]=PearCorrVarEst([ts(:,i),AdjCUR]',numel(AdjCUR),'taper','tukey',sqrt(T));
    
    [~,Stat]=xDF([ts(:,i),AdjCUR]',numel(AdjCUR),'truncate','adaptive','TVOn');
    
    %zstat(i) = Stat.z.rz(1,2);
    zMEfish(i) = Stat.z.rzf(1,2);
    
    %pstat(i) = Stat.p.r_Pval(1,2);
    pfish(i) = Stat.p.f_Pval(1,2);    

    %--AR1
    [~,~,~,bcf_ar1] = Bartlett46_fft([ts(:,i),AdjCUR]',numel(AdjCUR));
    z_ar1(i) = atanh(r(i)).*sqrt(numel(AdjCUR)./bcf_ar1(1,2)-3);

    %--C&R
    [~,~,~,bcf_CR] = CRBCF([ts(:,i),AdjCUR]',numel(AdjCUR));
    z_CR(i) = atanh(r(i)).*sqrt(numel(AdjCUR)./bcf_CR(1,2)-3);    
end

% atlasdir='/Users/sorooshafyouni/Home/HCP_Scripts/Yeo2011_17Networks_FSL_MNI152_2mm.nii.gz';
% PutMeOnYeo(atlasdir,znaive,'SPPOPGROWDEU')

%[~,~,pfish_adj] = fdr_bh(pfish,0.025);
%[~,~,p_adj] = fdr_bh(p,0.025);
%zMEfish_crtd = -z(p_adj);

p2z_crtd    = -z(FDR(p,alp));
p2z_ME_crtd = -z(FDR(pfish,alp));

sc_ME_rate = sum(fdr_bh(pfish))./114*100;

fpr_idx = sum(abs(znaive)>p2z & abs(zMEfish)<p2z)./114;
tpr_idx = find(pfish<alp);
tnr_idx = sum(abs(znaive)<p2z & abs(zMEfish)<p2z)./114;

fpr1_pidx = 85; 
tpr2_pidx = 37;

fpr2_pidx = 7;
tpr1_pidx = 66; % 87


r([fpr1_pidx,fpr2_pidx,tpr1_pidx,tpr2_pidx ]), 
disp('--')
znaive(fpr1_pidx)
p(fpr1_pidx)<alp, 
disp('--')
zMEfish(fpr1_pidx)
pfish(fpr1_pidx)<alp, 
%pstat(pidx)<alp

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%-----FPR1
%XCF
examplesfh = figure('position',[50,500,400,350]);
% subplot(4,2,1);
% hold on; grid on; box on;
% 
% [fpr1_xcf,fpr1_lags] = crosscorr(AdjCUR,ts(:,fpr1_pidx),T-1);
% plot(fpr1_lags,fpr1_xcf,'LineWidth',1.4,'Color',FP1col)
% 
% line([min(fpr1_lags) max(fpr1_lags)],[bnd bnd])
% line([min(fpr1_lags) max(fpr1_lags)],-[bnd bnd])
% 
% ylabel('Cross-correlation','fontsize',fs,'interpreter','latex')
% xlabel('Lags','fontsize',fs,'interpreter','latex')
% xlim([min(fpr1_lags) max(fpr1_lags)]/5)
% ylim([-1 1])

%------------------------------------------------------------------------
%ACF
% subplot(4,2,2);
% hold on; grid on; box on;
% plot(0:T-1,AC_RWV,'LineWidth',1.4,'Color',[1 0 0 .5])
% plot(0:T-1,AC_fft(ts(:,fpr1_pidx),T)','LineWidth',1.3,'Color',FP1col)
% 
% line([0 T-1], [bnd bnd])
% line([0 T-1],-[bnd bnd])
% 
% legend({'AC(UR)','AC(Node 85)'},'fontsize',fs,'location','northeast')
% ylabel('Autocorrelation','fontsize',fs,'interpreter','latex')
% xlabel('Lags','fontsize',fs,'interpreter','latex')
% xlim([0 T-1]/5)
% ylim([-1 1])

%-----FPR2
%------------------------------------------------------------------------
%XCF
subplot(2,2,1);
hold on; grid on; box on;

[fpr2_xcf,fpr2_lags] = crosscorr(AdjCUR,ts(:,fpr2_pidx),T-1);
plot(fpr2_lags,fpr2_xcf,'LineWidth',1.4,'Color',Blue)

line([min(fpr2_lags) max(fpr2_lags)],[bnd bnd])
line([min(fpr2_lags) max(fpr2_lags)],-[bnd bnd])

ylabel('Cross-correlation','fontsize',fs,'interpreter','latex')
xlabel('Lags','fontsize',fs,'interpreter','latex')
xlim([min(fpr2_lags) max(fpr2_lags)]/5)
ylim([-1 1])

%------------------------------------------------------------------------
%ACF
subplot(2,2,2);
hold on; grid on; box on;
plot(0:T-1,AC_RWV,'LineWidth',1.4,'Color',[0 0 0 .6])
plot(0:T-1,AC_fft(ts(:,fpr2_pidx),T)','LineWidth',1.3,'Color',Blue)

line([0 T-1],[bnd bnd])
line([0 T-1],-[bnd bnd])
%HCP118528
%HCP135932
legend({'L-PFCd','L-SoMotCent'},'fontsize',fs,'location','southwest')
ylabel('Autocorrelation','fontsize',fs,'interpreter','latex')
xlabel('Lags','fontsize',fs,'interpreter','latex')
%xlim([0 T-1]/5)
xlim([0 rightmxxlim])
ylim([-1 1])

%-----TPR1
%------------------------------------------------------------------------
%XCF
subplot(2,2,3);
hold on; grid on; box on;

[fpr1_xcf,fpr1_lags] = crosscorr(AdjCUR,ts(:,tpr1_pidx),T-1);
plot(fpr1_lags,fpr1_xcf,'LineWidth',1.4,'Color',Red)

line([min(fpr1_lags) max(fpr1_lags)],[bnd bnd])
line([min(fpr1_lags) max(fpr1_lags)],-[bnd bnd])

ylabel('Cross-correlation','fontsize',fs,'interpreter','latex')
xlabel('Lags','fontsize',fs,'interpreter','latex')
xlim([min(fpr1_lags) max(fpr1_lags)]/5)
ylim([-1 1])

%------------------------------------------------------------------------
%ACF
subplot(2,2,4);
hold on; grid on; box on;
plot(0:T-1,AC_RWV,'LineWidth',1.4,'Color',[0 0 0 .6])
plot(0:T-1,AC_fft(ts(:,tpr1_pidx),T)','LineWidth',1.3,'Color',Red)

line([0 T-1],[bnd bnd])
line([0 T-1],-[bnd bnd])
%HCP118528
%HCP135932
legend({'L-PFCd','R-Insula'},'fontsize',fs,'location','southwest')
ylabel('Autocorrelation','fontsize',fs,'interpreter','latex')
xlabel('Lags','fontsize',fs,'interpreter','latex')
%xlim([0 T-1]/5)
xlim([0 rightmxxlim])
ylim([-1 1])


%-----TPR2
%------------------------------------------------------------------------
%XCF
% subplot(4,2,5);
% hold on; grid on; box on;
% 
% [fpr2_xcf,fpr2_lags] = crosscorr(AdjCUR,ts(:,tpr2_pidx),T-1);
% plot(fpr2_lags,fpr2_xcf,'LineWidth',1.4,'Color',TP2col)
% 
% line([min(fpr2_lags) max(fpr2_lags)],[bnd bnd])
% line([min(fpr2_lags) max(fpr2_lags)],-[bnd bnd])
% 
% ylabel('Cross-correlation','fontsize',fs,'interpreter','latex')
% xlabel('Lags','fontsize',fs,'interpreter','latex')
% xlim([min(fpr2_lags) max(fpr2_lags)]/5)
% ylim([-1 1])

%------------------------------------------------------------------------
%ACF
% subplot(4,2,6);
% hold on; grid on; box on;
% plot(0:T-1,AC_RWV,'LineWidth',1.4,'Color',[1 0 0 .5])
% plot(0:T-1,AC_fft(ts(:,tpr2_pidx),T)','LineWidth',1.3,'Color',TP2col)
% 
% line([0 T-1],[bnd bnd])
% line([0 T-1],-[bnd bnd])
% 
% legend({'AC(UR)','AC(Node 37)'},'fontsize',fs,'location','northeast')
% ylabel('Autocorrelation','fontsize',fs,'interpreter','latex')
% xlabel('Lags','fontsize',fs,'interpreter','latex')
% xlim([0 T-1]/5)
% ylim([-1 1])
% 
set(examplesfh,'Color','w');
%------------------------------------------------------------------------
%------------------------------------------------------------------------

scatterfh=figure('position',[50,500,350,320]);
hold on; grid on; box on; %axis square
%title('HCP135392 v GPGR','FontSize',fs)

%patch([-xxlim -p2z -p2z -xxlim],[-p2z -p2z p2z p2z],[.5 .5 .5],'facealpha',0.2)
%patch(-[-xxlim -p2z -p2z -xxlim],[-p2z -p2z p2z p2z],[.5 .5 .5],'facealpha',0.2)

%dudu0=scatter(znaive,z_ar1,50,'marker','o','MarkerEdgeColor',Acol,'MarkerFaceColor',Acol,'MarkerFaceAlpha',0.30,'MarkerEdgeAlpha',0.30);
%dudu1=scatter(znaive,z_CR,50,'marker','d','MarkerEdgeColor',CRCol,'MarkerFaceColor',CRCol,'MarkerFaceAlpha',0.30,'MarkerEdgeAlpha',0.30); 
dudu2=scatter(znaive,zMEfish,50,'marker','s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.30,'MarkerEdgeAlpha',0.30); 
%dudu3=scatter(znaive,zMEfish_crtd,50,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.30,'MarkerEdgeAlpha',0.30); 
%legend([dudu0 dudu1 dudu2],{'Bartlett46','Clifford & Richardson','Monster Equation'})

xlabel('Naive Z-scores','FontSize',fs,'interpreter','latex')
ylabel('xDF Adjusted Z-scores','FontSize',fs,'interpreter','latex')

%On the y-axis 
lh_Nucorp =  line([p2z p2z],[-xxlim xxlim],'linewidth',lw,'linestyle',':');
             line(-[p2z p2z],[-xxlim xxlim],'linewidth',lw,'linestyle',':');
lh_Ncorp  =  line(-[p2z_crtd p2z_crtd],[-xxlim xxlim],'linewidth',lw);
             line([p2z_crtd p2z_crtd], [-xxlim xxlim] ,'linewidth',lw);
lh_Xuncorp = line([-xxlim xxlim],[p2z p2z],'color','r','linewidth',lw,'linestyle',':');
             line([-xxlim xxlim],-[p2z p2z],'color','r','linewidth',lw,'linestyle',':');
lh_Xcorp   = line([-xxlim xxlim],-[4 4],'color','r','linewidth',lw);
             line([-xxlim xxlim],[4 4],'color','r','linewidth',lw);

reflh = refline(1,0);
reflh.Color = 'k';
reflh.LineStyle = '-.';             
             
legend([lh_Nucorp lh_Ncorp lh_Xuncorp lh_Xcorp reflh],{'Uncorrected CV(Naive)','FDR CV(Naive)','Uncorrected CV(xDF)','FDR CV(xDF)','Ref. Line'},'fontsize',fs,'location','southeast')

ylim([-4.5 4.5])
xlim([-4.5 7])

%FP1
%scatter(znaive(fpr1_pidx),zMEfish(fpr1_pidx),60,'marker','s','MarkerEdgeColor','none','MarkerFaceColor',FP1col)
%scatter(znaive(fpr1_pidx),z_ar1(fpr1_pidx),60,'marker','o','MarkerEdgeColor','none','MarkerFaceColor',FP1col)
%scatter(znaive(fpr1_pidx),z_CR(fpr1_pidx),60,'marker','d','MarkerEdgeColor','none','MarkerFaceColor',FP1col)

%FP2
scatter(znaive(fpr2_pidx),zMEfish(fpr2_pidx),60,'marker','s','MarkerEdgeColor','none','MarkerFaceColor',Blue)
%scatter(znaive(fpr2_pidx),z_ar1(fpr2_pidx),60,'marker','o','MarkerEdgeColor','none','MarkerFaceColor',FP2col)
%scatter(znaive(fpr2_pidx),z_CR(fpr2_pidx),60,'marker','d','MarkerEdgeColor','none','MarkerFaceColor',FP2col)

%TP1
scatter(znaive(tpr1_pidx),zMEfish(tpr1_pidx),60,'marker','s','MarkerEdgeColor','none','MarkerFaceColor',Red)
%scatter(znaive(tpr1_pidx),z_ar1(tpr1_pidx),60,'marker','o','MarkerEdgeColor','none','MarkerFaceColor',TP1col)
%scatter(znaive(tpr1_pidx),z_CR(tpr1_pidx),60,'marker','d','MarkerEdgeColor','none','MarkerFaceColor',TP1col)

%TP2
%scatter(znaive(tpr2_pidx),zMEfish(tpr2_pidx),60,'marker','s','MarkerEdgeColor','none','MarkerFaceColor',TP2col)
%scatter(znaive(tpr2_pidx),z_ar1(tpr2_pidx),60,'marker','o','MarkerEdgeColor','none','MarkerFaceColor',TP2col)
%scatter(znaive(tpr2_pidx),z_CR(tpr2_pidx),60,'marker','d','MarkerEdgeColor','none','MarkerFaceColor',TP2col)

set(scatterfh,'Color','w');
%--------------------------------------------------------------------------
%BLANDATMAN----------------------------------------------------------------

%subplot(4,2,[5:8])
% figure; 
% hold on; grid on; box on; %axis square
% %title('HCP135392 v GPGR','FontSize',fs)
% ab_MEzfish  = abs(zMEfish);
% ab_zCR  = abs(z_CR);
% ab_zar1  = abs(z_ar1);
% ab_znaive = abs(znaive);
% 
% dudu0 = scatter((ab_zar1+ab_znaive)./2,ab_zar1-ab_znaive,'Marker','o','MarkerFaceColor',Acol,'MarkerFaceAlpha',.3,'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
% dudu1 = scatter((ab_zCR+ab_znaive)./2,ab_zCR-ab_znaive,'Marker','d','MarkerFaceColor',CRCol,'MarkerFaceAlpha',.3,'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
% dudu2 = scatter((ab_MEzfish+ab_znaive)./2,ab_MEzfish-ab_znaive,'Marker','s','MarkerFaceColor','k','MarkerFaceAlpha',.3,'MarkerEdgeColor','none','MarkerFaceAlpha',0.3);
% legend([dudu0 dudu1 dudu2],{'Bartlett35','Clifford&Richardson','xDF'})
% 
% %FP1
% %scatter((ab_zfish(fpr_idx)+ab_znaive(fpr_idx))./2,ab_zfish(fpr_idx)-ab_znaive(fpr_idx),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none');
% scatter((ab_MEzfish(fpr1_pidx)+ab_znaive(fpr1_pidx))./2,ab_MEzfish(fpr1_pidx)-ab_znaive(fpr1_pidx),60,'marker','s','MarkerEdgeColor','none','MarkerFaceColor',FP1col); 
% %scatter((ab_zar1(fpr1_pidx)+ab_znaive(fpr1_pidx))./2,ab_zar1(fpr1_pidx)-ab_znaive(fpr1_pidx),60,'marker','o','MarkerEdgeColor','none','MarkerFaceColor',FP1col);
% %scatter((ab_zCR(fpr1_pidx)+ab_znaive(fpr1_pidx))./2,ab_zCR(fpr1_pidx)-ab_znaive(fpr1_pidx),60,'marker','d','MarkerEdgeColor','none','MarkerFaceColor',FP1col); 
% 
% %FP2
% scatter((ab_MEzfish(fpr2_pidx)+ab_znaive(fpr2_pidx))./2,ab_MEzfish(fpr2_pidx)-ab_znaive(fpr2_pidx),60,'marker','s','MarkerEdgeColor','none','MarkerFaceColor',FP2col); 
% %scatter((ab_zar1(fpr2_pidx)+ab_znaive(fpr2_pidx))./2,ab_zar1(fpr2_pidx)-ab_znaive(fpr2_pidx),60,'marker','o','MarkerEdgeColor','none','MarkerFaceColor',FP2col); 
% %scatter((ab_zCR(fpr2_pidx)+ab_znaive(fpr2_pidx))./2,ab_zCR(fpr2_pidx)-ab_znaive(fpr2_pidx),60,'marker','d','MarkerEdgeColor','none','MarkerFaceColor',FP2col); 
% 
% %TP1
% scatter((ab_MEzfish(tpr1_pidx)+ab_znaive(tpr1_pidx))./2,ab_MEzfish(tpr1_pidx)-ab_znaive(tpr1_pidx),60,'marker','s','MarkerEdgeColor','none','MarkerFaceColor',TP1col); 
% %scatter((ab_zar1(tpr1_pidx)+ab_znaive(tpr1_pidx))./2,ab_zar1(tpr1_pidx)-ab_znaive(tpr1_pidx),60,'marker','o','MarkerEdgeColor','none','MarkerFaceColor',TP1col); 
% %scatter((ab_zCR(tpr1_pidx)+ab_znaive(tpr1_pidx))./2,ab_zCR(tpr1_pidx)-ab_znaive(tpr1_pidx),60,'marker','d','MarkerEdgeColor','none','MarkerFaceColor',TP1col); 
% 
% %TP2
% scatter((ab_MEzfish(tpr2_pidx)+ab_znaive(tpr2_pidx))./2,ab_MEzfish(tpr2_pidx)-ab_znaive(tpr2_pidx),60,'marker','s','MarkerEdgeColor','none','MarkerFaceColor',TP2col); 
% %scatter((ab_zar1(tpr2_pidx)+ab_znaive(tpr2_pidx))./2,ab_zar1(tpr2_pidx)-ab_znaive(tpr2_pidx),60,'marker','o','MarkerEdgeColor','none','MarkerFaceColor',TP2col); 
% %scatter((ab_zCR(tpr2_pidx)+ab_znaive(tpr2_pidx))./2,ab_zCR(tpr2_pidx)-ab_znaive(tpr2_pidx),60,'marker','d','MarkerEdgeColor','none','MarkerFaceColor',TP2col); 
% 
% 
% line([0 10],[0 0],'linewidth',1.2,'color','k','linestyle','-.')
% xlim([0 4])
% ylabel('|xDF z-score| - |Fisher z-score|','FontSize',fs,'interpreter','latex')
% xlabel('<xDF z-score,Fisher z-score>','FontSize',fs,'interpreter','latex')


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%figure; hold on; box on; 
%scatter(pfish_adj,p_adj)

% figure('position',[50,500,350,350]);
% hold on; grid on; box on; axis square
% scatter(ac_strength,ab_zfish-ab_znaive,90,'marker','.','MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
% corr(ac_strength,(ab_zfish-ab_znaive)')
% refline
% xlabel('Autocorrelation Strength','FontSize',fs)
% ylabel('Changes in z-scores','FontSize',fs)

%------------------------------------------------------------------------
%------------------------------------------------------------------------

ts0 = ts(:,fpr2_pidx)-mean(ts(:,fpr2_pidx));
ts1 = ts(:,tpr1_pidx)-mean(ts(:,tpr1_pidx));

ts0 = ts0./std(ts0);
ts1 = ts1./std(ts1);


BOLDResfh = figure('position',[50,500,1100,200]);
hold on; grid on; box on; 
title('HCP135932')
plot(timspan,ts1,'color',[Red 0.5],'linewidth',lw)
plot(timspan,ts0,'color',[Blue 0.9],'linewidth',lw)
set(BOLDResfh,'color','w')
xlim([0 scanses])
legend({'Right Insula (R-Insula)','Left SomatoMotor Central (L-SomMot Cent)'},'fontsize',fs)
ylabel('Standardised BOLD Response','fontsize',fs,'interpreter','latex')
xlabel('Seconds','fontsize',fs,'interpreter','latex')
%------------------------------
% 
 
export_fig(examplesfh,'Figs/Examples.pdf') 
export_fig(ph0,'Figs/ph0.pdf') 
export_fig(scatterfh,'Figs/FisherPlot.pdf')
export_fig(BOLDResfh,'Figs/BOLDRes.pdf')