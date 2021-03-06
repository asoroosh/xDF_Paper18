%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7:
% Plots changes on the global measures of Statistically thresholded 
% (FDR-corrected) networks due to xDF vs. BH corrections on Z-scores. 
% The results were shown using Bland-Altman plots. 
%
% REQUIREMENTS:
% Data: Graph measures estimated for HCP 100 unrelated package
% 
% Code: 
%
% Soroosh Afyouni, University of Oxford, 2019, 
% srafyouni@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
%_____________________________________________________
load('/Users/sorooshafyouni/Home/BCF/BCFAnal/FCAnal/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;
SubList(3) = [];

VarList={'Den','GE'}; %,'C'};
VarList2Show={'\textbf{FDR Density}','\textbf{GE}'};

Col=get(groot,'defaultAxesColorOrder');
%------

AtlasList={'Power','Yeo','Gordon','ICA200'};%,'Gordon'};
fdMList={'Naive','MEs','CR'};
fdMList_labels={'Naive','xDF','BH'};

YLims = {[0 0.45],[0 8]};

PP = 'FPP';

thrM_L={'fdr'};

load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Yeo17_Networks_FuncBlkAss.mat'     ,'YeoZ','YeoU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Power2011_Networks_FuncBlkAss.mat' ,'PwrZ','PwrU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Gordon_13Networks_FuncBlkAss.mat'  ,'GordonZ','GordonU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/AAL2_Networks_FuncBlkAss.mat'      ,'AALZ','AALU')

Z={PwrZ;YeoZ;GordonZ};
U={PwrU;YeoU;GordonU};

cnt_g=1; spc=1;

fh = figure('position',[50,500,700,260]);
hold on; box on;

cnt_v=1;
for TarVar = VarList
    cnt_a=1;
    for Atlas = AtlasList
            disp(['A: ' Atlas{1}])
            for s_cnt = 1:numel(SubList)
                GMDir=['/Users/sorooshafyouni/Home/BCF/BCFAnal/FCAnal/R_GM/Stat/' Atlas{1} '/' PP '/'];

                Naive = load([GMDir 'GM_' Atlas{1} '_' PP '_' fdMList{1} '_' SubList{s_cnt} '_TVOff.mat']);
                ME   = load([GMDir 'GM_' Atlas{1} '_' PP '_' fdMList{2} '_' SubList{s_cnt}  '_TVOff.mat']);
                CR   = load([GMDir 'GM_' Atlas{1} '_' PP '_' fdMList{3} '_' SubList{s_cnt}  '_TVOff.mat']);
                
                ndf(s_cnt,1)   = eval(['Naive.w.' TarVar{1}]);
                hmedf(s_cnt,1) = eval(['ME.w.'   TarVar{1}]);
                hcrdf(s_cnt,1) = eval(['CR.w.'   TarVar{1}]);
                %__________________________________________________________________
            end
            
            %__________________________________________________________________           
            sph = subplot(length(VarList),length(AtlasList),spc);
            hold on; box on; 
            title([VarList2Show{cnt_v} '(' Atlas{1} ')'],'fontsize',12,'interpreter','latex')
            %__________________________________________________________________
            hndls = ScatterBoxPlots([ndf,hmedf,hcrdf],'subplot',sph,'color',Col(1:3,:),'pointsize',15,'MedLineColor','k');
            
            sph.XTickLabel=fdMList_labels;
            sph.FontSize = 11;
            sph.XTickLabelRotation=45;
            
            ylim(YLims{cnt_v});
            
            Htest(cnt_v,cnt_a) = ttest2(hmedf,hcrdf,'vartype','unequal');
            
            clear ndf hndf Naive Adj
            cnt_a = cnt_a+1;
            spc   = spc+1;
    end
    
    cnt_v = cnt_v+1;
end

set(gcf,'Color','w');

export_fig(fh,'Global_FDRStatDen.pdf');