clear

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/DVARS'))
addpath(genpath('/Users/sorooshafyouni/Home/GitClone/HetBiv'))

fs = 12;

load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/S/HCP_100Unrel_SubList.mat'])
SubList = HCP_10Unrel_SubList;
SubList(3) = [];

ROIList = {'LRPCC','LRSOMMOT','SVNPFCd'};
ROILabels = {'LH-PCC','LH-SomMot','LH-PFCd'};
ROIIdx = {[47 105],[6 63],[24 82]};

T = 1200; 
fh0 = figure('position',[50,500,410,400]);

s_cnt = 1;
for s=SubList
    r_cnt = 1; 
    for r = ROIList
        load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/ROIsVSCorrLeng/R/HCP100UR/' r{1} 'ACLs/HCP_' s{1} '_LR_' r{1} '_ACL.mat'],'ROIACL')
        
        var_Voxel_ACL(:,r_cnt,s_cnt) = [var(ROIACL{1}) var(ROIACL{2})];
        Voxel_ACL(:,r_cnt,s_cnt) = [mean(ROIACL{1}) mean(ROIACL{2})]; 

        load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/FC/100HCPTimeSeries/Yeo/HCP_FPP_' s{1} '_OnlyMTS.mat'])
        
        mts = mts - mean(mts);
        mts = GSRme(mts,T);
        
        ROI_ACL(:,r_cnt,s_cnt) = sum(AC_fft(mts(:,ROIIdx{r_cnt}),T).^2,2)';
        
        clear AC
        r_cnt = r_cnt + 1; 
    end
    s_cnt = s_cnt + 1;
end

figure(fh0)
sp0 = subplot(2,2,1);
YYLIM = 12.5; 
hold on; grid on; box on; 
title('Voxel Time Series','fontsize',fs,'interpreter','latex')
ScatterBoxPlots(squeeze(Voxel_ACL(1,:,:))','subplot',sp0,'pointsize',15);
ylim([1 YYLIM])
sp0.YTick = 0:2:YYLIM;
sp0.XTick = 1:3;
sp0.XTickLabel = ROILabels;
sp0.XTickLabelRotation = 45; 
ylabel('Autocorrelation Length','fontsize',fs,'interpreter','latex')

sp1 = subplot(2,2,2);
%YYLIM = 15;
hold on; grid on; box on; 
title('Averaged Time Series','fontsize',fs,'interpreter','latex')
ScatterBoxPlots(squeeze(ROI_ACL(1,:,:))','subplot',sp1,'pointsize',15);
ylim([1 YYLIM])
sp1.YTick = 0:2:YYLIM;
sp1.XTick = 1:3;
sp1.XTickLabel = ROILabels;
sp1.XTickLabelRotation = 45; 

ylabel('Autocorrelation Length','fontsize',fs,'interpreter','latex')

set(fh0,'Color','w');
