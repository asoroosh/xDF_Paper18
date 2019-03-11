%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FisherPlot.m
%
% Generates Figure 6. 
% This script compares the differences in FC strength after autocorrelation
% corrections. 
% 
%%%%CODE DEPENDENCIES:
%   xDF.m
%   CRBCF.m
%   fdr_bh.m 
%   CostEff_bin.m 
%   
%
%%%%DATA DEPENDENCIES:
%   You need to have the time series of subject 118528 and/or 135932 of the
%   HCP cohort parcellated with Yeo2001 atlas. 
%   
%   Data should be openly available via https://db.humanconnectome.org
%
%%%%HOW TO RUN
%  Change the path directories for addpath and load lines and click Run!
%
%   Please contact srafyouni@gmail.com for any questions/bugs. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soroosh Afyouni, University of Oxford, 2018

clear

addpath /Users/sorooshafyouni/Home/matlab/Ext/StatDens
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/xDF'))
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/DVARS'))

fs=12;
lw=1.3;

xxlim=50;

z = @(p) -sqrt(2) .* erfcinv(p.*2);
alp=0.05;

z2p = -z(alp/2);

Col=get(groot,'defaultAxesColorOrder');
FP1col  = Col(1,:); % Blue
FP2col  = Col(3,:); % Yellow
TP1col  = Col(2,:); % Red (/orange!)
TP2col  = Col(4,:); % Purple
Acol    = Col(5,:); % Green

%==========================================================================
% SubID = 135932;
% Wnodes = [37 87];
% Anodes = [7 85];
% Fnodes = [88 23];

SubID = 118528;
Wnodes = [94 37];
Anodes = [25 13];
Fnodes = [104 103];

load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_' num2str(SubID) '_OnlyMTS.mat'],'mts')
%==========================================================================

%----For Cost Efficient Density Detection:
dnsrng = 0.05:0.025:0.50;
OnlyPositives = 1;
Binarise = 0; 
%----

T  = size(mts,1);
nn = size(mts,2);

%----GSR ME
mts = GSRme(mts,T); %Global Signal Regression

Idx = find(triu(ones(nn),1));
ne  = numel(Idx);

%prepare the Naive -------------------------------------------------
[rmat,P_naive]=corr(mts);

Z_naive = atanh(rmat).*sqrt(T-3);
Z_naive = real(Z_naive);
Z_naive_triu = Z_naive(Idx);
Z_naive(1:nn+1:end) = 0;

P_naive_triu = P_naive(Idx);
[~,CP_naive,FDRAdjP_naive_triu] = fdr_bh(P_naive_triu);
Z_cv_Naive = -z(CP_naive./2);

[~,CE_den_naive,~] = CostEff_bin(Z_naive,dnsrng,OnlyPositives,Binarise);
Z_ce_Naive = prctile(Z_naive_triu,(1-CE_den_naive)*100);
%prepare the monster-------------------------------------------------
disp('Var is estimated...')
[~,Stat] = xDF(mts,T,'truncate','adaptive','TVOff','verbose');
Z_ME = Stat.z; 
P_ME = Stat.p; 
Z_ME         = real(Z_ME);
Z_ME_triu    = Z_ME(Idx);
P_ME_triu    = P_ME(Idx);

[~,CP_ME,FDRAdjP_ME_triu]   = fdr_bh(P_ME_triu);
Z_cv_ME             = -z(CP_ME./2);

[~,CE_den_ME,~] = CostEff_bin(Z_ME,dnsrng,OnlyPositives,Binarise);
Z_ce_ME = prctile(Z_ME_triu,(1-CE_den_ME)*100);

%prepare the CR------------------------------------------------------
[~,Z_CR,P_CR] = CRBCF(mts,T);
%z_CR = atanh(rmat).*sqrt((T./CR)-3);
%z_CR = real(z_CR);
Z_CR_triu   = Z_CR(Idx);
P_CR_triu   = P_CR(Idx);

[~,CP_CR,FDRAdjP_CR_triu]   = fdr_bh(P_CR_triu);
Z_cv_CR             = -z(CP_CR./2);

[~,CE_den_CR,~] = CostEff_bin(Z_CR,dnsrng,OnlyPositives,Binarise);
Z_ce_CR = prctile(Z_CR_triu,(1-CE_den_CR)*100);

%-----------------------------------------------------------------------------
%FISHER PLOTS-----------------------------------------------------------------
%-----------------------------------------------------------------------------

fisherplots = figure('position',[50,500,1350,330]);
figure(fisherplots)
subplot(1,3,1)
hold on; grid off; box on; axis square
%patch([-xxlim -z2p -z2p -xxlim],[-z2p -z2p z2p z2p],[.5 .5 .5],'facealpha',0.2)
%patch(-[-xxlim -z2p -z2p -xxlim],[-z2p -z2p z2p z2p],[.5 .5 .5],'facealpha',0.2)

scatter(Z_naive_triu,Z_ME_triu,50,'marker','.','MarkerEdgeColor',TP1col,'MarkerFaceColor',TP1col,'MarkerEdgeAlpha',0.6);

FDRME_Idx = (abs(Z_naive_triu)>Z_cv_Naive  & abs(Z_ME_triu)<Z_cv_ME); %|| (Z_naive_triu>-Z_cv_Naive  && Z_ME_triu<-Z_cv_ME);
CEME_Idx =  (abs(Z_naive_triu)>Z_ce_Naive  & abs(Z_ME_triu)<Z_ce_ME); %|| (Z_naive_triu>-Z_ce_Naive  && Z_ME_triu<-Z_ce_ME); 

FPRME_dots = scatter(Z_naive_triu(FDRME_Idx),Z_ME_triu(FDRME_Idx),50,'marker','.','MarkerEdgeColor',Acol,'MarkerFaceColor',Acol,'MarkerEdgeAlpha',1);
CEME_dots = scatter(Z_naive_triu(CEME_Idx),Z_ME_triu(CEME_Idx),50,'marker','.','MarkerEdgeColor',FP1col,'MarkerFaceColor',FP1col,'MarkerEdgeAlpha',1);

Wnodes_dots = scatter(Z_naive(Wnodes(1),Wnodes(2)),Z_ME(Wnodes(1),Wnodes(2)),40,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
Anodes_dots = scatter(Z_naive(Anodes(1),Anodes(2)),Z_ME(Anodes(1),Anodes(2)),40,'marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');

fpr_ucor_ME = sum(abs(Z_naive_triu)>z2p        & abs(Z_ME_triu)<z2p)./ne;
fpr_fdr_ME  = sum(FDRME_Idx)./ne;
fpr_ce_ME  = sum(CEME_Idx)./ne;

fpr_ucor_CR = sum(abs(Z_naive_triu)>z2p         & abs(Z_CR_triu)<z2p)./ne;
fpr_fdr_CR  = sum(abs(Z_naive_triu)>Z_cv_Naive  & abs(Z_CR_triu)<Z_cv_CR)./ne;
fpr_ce_CR   = sum(abs(Z_naive_triu)>Z_ce_Naive  & abs(Z_CR_triu)<Z_ce_CR)./ne;

xlabel(['Naive Z-scores'],'FontSize',fs,'Interpreter','latex')
ylabel(['xDF Adjusted Z-scores'],'FontSize',fs,'Interpreter','latex')

%Uncorrected critical values
%line([z2p z2p],[-xxlim xxlim],'color' ,'b','linestyle',':','linewidth',1.2)
%line(-[z2p z2p],[-xxlim xxlim],'color','b','linestyle',':','linewidth',1.2)
%line([-xxlim xxlim],[z2p z2p],'color' ,'r','linestyle',':','linewidth',1.2)
%line([-xxlim xxlim],-[z2p z2p],'color','r','linestyle',':','linewidth',1.2)

%FDR-corrected critical values
Line_fdrOnX_Naive = line([Z_cv_Naive Z_cv_Naive],[-xxlim xxlim],'color' ,[0 0 1 0.7],'linestyle','-.','linewidth',1.2); %on naive; x-axis
                    line(-[Z_cv_Naive Z_cv_Naive],[-xxlim xxlim],'color',[0 0 1 0.7],'linestyle','-.','linewidth',1.2)
Line_fdrOnY_ME =    line([-xxlim xxlim],[Z_cv_ME Z_cv_ME],'color' ,[1 0 0 0.7],'linestyle','-.','linewidth',1.2); %on Monster Eq.; x-axis
                    line([-xxlim xxlim],-[Z_cv_ME Z_cv_ME],'color',[1 0 0 0.7],'linestyle','-.','linewidth',1.2)    

%Cost Efficient critical values
Line_CEOnX_Naive =  line([Z_ce_Naive Z_ce_Naive],[-xxlim xxlim],'color' ,[0 0 1 0.7],'linestyle','-','linewidth',1.2); %on naive; x-axis
                    line(-[Z_ce_Naive Z_ce_Naive],[-xxlim xxlim],'color',[0 0 1 0.7],'linestyle','-','linewidth',1.2)
Line_CEOnY_ME =     line([-xxlim xxlim],[Z_ce_ME Z_ce_ME],'color' ,[1 0 0 0.7],'linestyle','-','linewidth',1.2); %on Monster Eq.; x-axis
                    line([-xxlim xxlim],-[Z_ce_ME Z_ce_ME],'color',[1 0 0 0.7],'linestyle','-','linewidth',1.2)        

%blah blah blah
%ylim([-xxlim xxlim])
%xlim([-xxlim xxlim])

ylim([-10 10])
xlim([-10 10])

reflh = refline(1,0);
reflh.LineWidth = 1.4;
reflh.LineStyle = '-.';
reflh.Color = 'k';

legend([Line_fdrOnX_Naive Line_fdrOnY_ME Line_CEOnX_Naive Line_CEOnY_ME reflh FPRME_dots CEME_dots],...
    {'FDR CV(Naive)','FDR CV(xDF)','CE(Naive)','CE(xDF)','Ref. Line','FP(FDR)','FP(CE)'},'fontsize',fs,'location','northwest')

%--------------------------------------------------------------------------
subplot(1,3,2)
hold on; grid off; box on; axis square; grid off;
%patch([-xxlim -z2p -z2p -xxlim],[-z2p -z2p z2p z2p],[.5 .5 .5],'facealpha',0.2)
%patch(-[-xxlim -z2p -z2p -xxlim],[-z2p -z2p z2p z2p],[.5 .5 .5],'facealpha',0.2)

scatter(Z_CR_triu,Z_ME_triu,50,'marker','.','MarkerEdgeColor',TP2col,'MarkerFaceColor',TP2col,'MarkerEdgeAlpha',0.9); 

Wnodes_dots = scatter(Z_CR(Wnodes(1),Wnodes(2)),Z_ME(Wnodes(1),Wnodes(2)),40,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');
Anodes_dots = scatter(Z_CR(Anodes(1),Anodes(2)),Z_ME(Anodes(1),Anodes(2)),40,'marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');
Cnodes_dots = scatter(Z_CR(Fnodes(1),Fnodes(2)),Z_ME(Fnodes(1),Fnodes(2)),40,'marker','o','MarkerEdgeColor','g','MarkerFaceColor','g');

FDRMECR_Idx = abs(Z_CR_triu)<Z_cv_CR  & abs(Z_ME_triu)>Z_cv_ME;
UCMECR_Idx  = abs(Z_CR_triu)>z2p      & abs(Z_ME_triu)<z2p;
CEMECR_Idx  = abs(Z_CR_triu)<Z_ce_CR  & abs(Z_ME_triu)>Z_ce_ME;

fpr_ucor_MECR   = sum(UCMECR_Idx)./ne;
fpr_fdr_MECR    = sum(FDRMECR_Idx)./ne; %thos edges where detected back...
fpr_ce_MECR     = sum(CEMECR_Idx)./ne;

xlabel(['BH adjusted Z-scores'],'FontSize',fs,'Interpreter','latex')
ylabel(['xDF Adjusted Z-scores'],'FontSize',fs,'Interpreter','latex')


%FDR-corrected critical values
FDRCR_lh =  line([Z_cv_CR Z_cv_CR],[-xxlim xxlim],'color' ,[0 0 1 0.7],'linestyle','-.','linewidth',1.2); %on naive; x-axis
            line(-[Z_cv_CR Z_cv_CR],[-xxlim xxlim],'color',[0 0 1 0.7],'linestyle','-.','linewidth',1.2);
FDRME_lh =  line([-xxlim xxlim],[Z_cv_ME Z_cv_ME],'color' ,[1 0 0 0.7],'linestyle','-.','linewidth',1.2); %on Monster Eq.; x-axis
            line([-xxlim xxlim],-[Z_cv_ME Z_cv_ME],'color',[1 0 0 0.7],'linestyle','-.','linewidth',1.2);    

%Cost Efficient critical values
CECR_lh =   line([Z_ce_CR Z_ce_CR],[-xxlim xxlim],'color' ,[0 0 1 0.7],'linestyle','-','linewidth',1.2); %on naive; x-axis
            line(-[Z_ce_CR Z_ce_CR],[-xxlim xxlim],'color',[0 0 1 0.7],'linestyle','-','linewidth',1.2);
CRME_lh =   line([-xxlim xxlim],[Z_ce_ME Z_ce_ME],'color' ,[1 0 0 0.7],'linestyle','-','linewidth',1.2); %on Monster Eq.; x-axis
            line([-xxlim xxlim],-[Z_ce_ME Z_ce_ME],'color',[1 0 0 0.7],'linestyle','-','linewidth',1.2);      

ylim([-20 20])
xlim([-20 20])
reflh = refline(1,0);
reflh.LineWidth = 1.4;
reflh.LineStyle = '-.';
reflh.Color = 'k';    

legend([FDRCR_lh FDRME_lh CECR_lh CRME_lh reflh],...
    {'FDR CV(BH)','FDR CV(xDF)','CE(BH)','CE(xDF)','Ref. Line'},'fontsize',fs,'location','northwest')

%-----------------------------------------------------------------------------
%PVALS------------------------------------------------------------------------
%-----------------------------------------------------------------------------

%figure(pvalsfh);
sp0 = subplot(1,3,3); 
hold on; box on; grid on; axis square; 

% %scatter(P_naive_triu,pMEfish_triu,60,'marker','.','MarkerEdgeColor',FP2col,'MarkerFaceColor',FP2col,'MarkerEdgeAlpha',0.5)
% %scatter(P_CR_triu,pMEfish_triu   ,60,'marker','.','MarkerEdgeColor',TP2col,'MarkerFaceColor',TP2col,'MarkerEdgeAlpha',0.1)
% 
% %set(sp0,'defaultAxesColorOrder',[FP2col ; TP2col])
% 
% yyaxis left
% scatter(FDRAdjP_ME_triu,FDRAdjP_naive_triu,60,'marker','.','MarkerEdgeColor',TP1col,'MarkerFaceColor',TP1col,'MarkerEdgeAlpha',.9)
% ylabel('Naive Adjusted p-values (FDR-corrected)','Interpreter','latex','FontSize',fs)
% set(gca,'YColor',TP1col)
% yyaxis right
% scatter(FDRAdjP_ME_triu,FDRAdjP_CR_triu   ,60,'marker','.','MarkerEdgeColor',TP2col,'MarkerFaceColor',TP2col,'MarkerEdgeAlpha',.9)
% ylabel('BH Adjusted p-values (FDR-corrected)','Interpreter','latex','FontSize',fs)
% xlabel('xDF Adjusted p-values (FDR-corrected)','Interpreter','latex','FontSize',fs)
% set(gca,'YColor',TP2col)
% %legend({'p-values (xNaive) ',''})
% 
% reflh = refline(1,0);
% reflh.LineWidth = 1.4;
% reflh.Color = 'k';  
% reflh.LineStyle = '-.';
% %ylabel('ME p-values','fontsize',fs)
% %xlabel('p-values','fontsize',fs)
% 
% set(fisherplots,'Color','w');
% 
% % export_fig(zstatfh,'Figs/zstat_FC.pdf')
% % 
% % set(fisherplots,'Color','w');
% % export_fig(fisherplots,'Figs/FC_FisherPlots.pdf')
% % 
% % set(pvalsfh,'Color','w');

BonPvals = -log(0.05./numel(P_ME_triu));

plot(-log(sort(P_ME_triu,'descend')),'linewidth',1.4,'color',Acol); 
plot(-log(sort(P_CR_triu,'descend')),'linewidth',1.4,'color',FP1col); 
plot(-log(sort(P_naive_triu,'descend')),'linewidth',1.4,'color',TP1col);
line([0 numel(P_naive_triu)],[BonPvals BonPvals],'linewidth',1.4,'color',[0 0 0 0.5],'linestyle','-.')
ylabel('-log(p-values)','Interpreter','latex','FontSize',fs)
xlabel('p-value Order','Interpreter','latex','FontSize',fs)

legend({'xDF','BH','Naive','Bonferroni -log(p-value)'},'fontsize',fs,'location','southeast')

ylim([0 15])
xlim([0 numel(P_naive_triu)])

set(fisherplots,'Color','w');

%-----------------------------------------------------------------------------
%Autocorrelations-------------------------------------------------------------
%-----------------------------------------------------------------------------
[ac,bnd] = AC_fft(mts,T);
bnd = bnd(2);


fh = figure('position',[50,500,600,200]);
subplot(1,2,1)
grid on; box on; 
plot(ac([88 23],1:100)','marker','.','linewidth',lw);
legend({'AC(Node88)','AC(Node23)'},'fontsize',fs)
ylabel('Autocorrelation','fontsize',fs,'interpreter','latex')
xlabel('Lags','fontsize',fs,'interpreter','latex')

line([0 100],[bnd bnd])
line([0 100],-[bnd bnd])
xlim([0 100])

subplot(1,2,2)
grid on; box on; 

[fpr1_xcf,fpr1_lags] = crosscorr(mts(:,88),mts(:,23),T-1);
plot(fpr1_lags,fpr1_xcf,'marker','.','linewidth',lw);
legend({'XC(Node88,Node23)'},'fontsize',fs)
ylabel('Cross-correlation','fontsize',fs,'interpreter','latex')
xlabel('Lags','fontsize',fs,'interpreter','latex')

line([0 100],[bnd bnd])
line([0 100],-[bnd bnd])
xlim([0 100])

set(gcf,'Color','w');

%export_fig(fisherplots,['Figs/FisherPlotsOfHCPSub_' num2str(SubID) '_Amended.pdf'])
%export_fig(fh,['Figs/AutoCrossCrorr_Fnodes_' num2str(SubID) '.pdf'])

disp('----')
disp(['FPR FDR ME: ', num2str(fpr_fdr_ME*100)])
disp(['FPR FDR BH: ', num2str(fpr_fdr_CR*100)])
disp('----')
disp(['FPR CE ME: ', num2str(fpr_ce_ME*100)])
disp(['FPR CE BH: ', num2str(fpr_ce_CR*100)])



