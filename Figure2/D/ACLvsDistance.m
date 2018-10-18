clear

load(['/Volumes/Extra/HCP_100UR_FC/Power/FPP/HCP_FPP_100307_ROIsTS.mat'],'ROI')
N = numel(ROI.ts);

Col = get(groot,'defaultAxesColorOrder');

%Deep_loc = [44 45 39];
Deep_loc = [45 57 38];
%Deep_loc = mean(ROI.XYZ{222});

load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Power2011_Networks_FuncBlkAss.mat' ,'PwrZ','PwrU')
%TemporalLobes = ismember(PwrZ,[247,248,81,81,249,250,128,129]);
CortSubCortIdx = ones(1,N);
CortSubCortIdx(find(PwrZ==10 | PwrZ==11 )) = 0; 

load('/Users/sorooshafyouni/Home/DVARS/fMRIDiag/HCP/GordonParcellations/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;

SubList(3) = [];

for i = 1:264
    ROI_loc = mean(ROI.XYZ{i});
    ROI_sz(i) = size(ROI.XYZ{i},1);
    ROIdist(i) = sqrt(sum((ROI_loc-Deep_loc).^2));
end

load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/ROIsVSCorrLeng/R/HCP100UR/Power_GSR_FPP_CorrLeng_ROIWise.mat'])

% T = 1200;
% for i = 1:99
%     load(['/Volumes/Extra/HCP_100UR_FC/Power/FPP/HCP_FPP_' SubList{i} '_ROIsTS.mat'],'mts');
%     
%     AC_tmp = AC_fft(mts,T);
%     AC(:,i) = sum(abs(AC_tmp(:,2:15)),2);
%     
% end

mACL = mean(ACL,2);

SNR_AC = corr(mACL,ROIdist');

%[xl,yl]=FitMeLine(mean(ACL,2),ROIdist',2);

fh = figure('position',[50,500,400,350]); 
hold on; box on; grid on;
% for i = 1:99
%    ACLDist(i) = corr(ACL5(:,i),ROIdist');
%    scatter(ACL5(:,i),sqrt(ROIdist),5,[.5 .5 .5]) 
% end

gscatter(ROIdist,mean(ACL,2),~CortSubCortIdx,Col([3 1],:));

text(10,5,['R= ' num2str(round(SNR_AC,2))],'FontSize',12)
text(10,4.5,['R^2= ' num2str(round(SNR_AC.^2,2))],'FontSize',12)

legend({'Cortical','Subcortical \& Cerebellar'},'fontsize',12,'Interpreter','latex')
%scatter(ROIdist,mean(ACL,2),'MarkerFaceColor',Col(3,:),'MarkerEdgeColor',Col(3,:))

xlim([0 50]); ylim([1 6.5])
%refline(1,0)

ylabel('Autocorrelation Length')
xlabel(['Euclidean Distance to [' num2str(Deep_loc) '] (voxels)'])

%plot(xl,yl)

set(fh,'Color','w');
export_fig(fh,'ACLvDistance.pdf')
