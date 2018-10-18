clear

fs = 12;
lw = 1.3;

%----------------

load('/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_118528_OnlyMTS.mat','mts')
ac = AC_fft(mts,1200);
ac(:,[1 end]) = [];

acfh = figure('position',[50,500,600,150]); 
hold on; box on; axis tight
imagesc(abs(ac(:,1:200)),[0 1])
ylabel('Nodes','fontsize',fs,'interpreter','latex')
xlabel('Lags','fontsize',fs,'interpreter','latex')
colormap gray
set(acfh,'color','w')

export_fig(acfh,'Figs/Sim_118528.pdf')

%----------------

load('/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_135932_OnlyMTS.mat','mts')
ac = AC_fft(mts,1200);
ac(:,[1 end]) = [];

acfh = figure('position',[50,500,600,150]); 
hold on; box on; axis tight
imagesc(abs(ac(:,1:200)),[0 1])
ylabel('Nodes','fontsize',fs,'interpreter','latex')
xlabel('Lags','fontsize',fs,'interpreter','latex')
colormap gray
set(acfh,'color','w')

export_fig(acfh,'Figs/Sim_135932.pdf')
%----------------

acfh = figure('position',[50,500,200,150]); 
hold on; box on; grid on; 
plot(ac(35,1:20),'linewidth',lw,'marker','.','color','k')
xlabel('Lags','fontsize',fs,'interpreter','latex')
ylabel('Autocorrelation Coefficient','fontsize',fs,'interpreter','latex')
set(acfh,'color','w')

export_fig(acfh,'Figs/SimAC.pdf')
%----------------

load('/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_118528_OnlyMTS.mat','mts')
mat = corr(mts);
mat = threshold_proportional(mat,0.15);
mat(mat>0) = 1;

acfh = figure('position',[50,500,200,200]);
hold on; axis square; axis tight; 
imagesc(mat)
ylabel('Nodes','fontsize',fs,'interpreter','latex')
xlabel('Nodes','fontsize',fs,'interpreter','latex')
colormap gray
set(acfh,'color','w')

export_fig(acfh,'Figs/SimMat.pdf')
