
clear 
%_____________________________________________________
load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FCAnal/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;
SubList(3) = [];

VarList={'fdr','CE'}; %,'C'};
VarList2Show={'\textbf{wCE}','\textbf{FDR}'};

Col=get(groot,'defaultAxesColorOrder');
%------

fs = 12;

AtlasList={'Power','Yeo','Gordon'};%,'Gordon'};
%fdMList={'Naive','MEs','CR'};
fdMList_labels={'Naive','xDF','BH'};

YLims = {[3 15],[3 15]};

PP = 'FPP';

thrM_L={'fdr'};

load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Yeo17_Networks_FuncBlkAss.mat'     ,'YeoZ','YeoU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Power2011_Networks_FuncBlkAss.mat' ,'PwrZ','PwrU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Gordon_13Networks_FuncBlkAss.mat'  ,'GordonZ','GordonU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/AAL2_Networks_FuncBlkAss.mat'      ,'AALZ','AALU')

Z={PwrZ;YeoZ;GordonZ};
U={PwrU;YeoU;GordonU};

cnt_g=1; spc=1;
fh = figure('position',[50,500,600,300]);
hold on; box on;    
cnt_v=1;
for TarVar = VarList
    cnt_a=1;
    for Atlas = AtlasList
            disp(['A: ' Atlas{1}])
            for s_cnt = 1:numel(SubList)
                GMDir=['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FCAnal/R_mats/' Atlas{1} '/' PP '/'];

                load([GMDir 'Mats_' Atlas{1} '_' PP '_'  SubList{s_cnt} '.mat']);
                
                ndf(s_cnt,1)   = mean(nonzeros(eval(['twMats.Naive_' TarVar{1}])));
                hmedf(s_cnt,1) = mean(nonzeros(eval(['twMats.ME_' TarVar{1}])));
                hcrdf(s_cnt,1) = mean(nonzeros(eval(['twMats.CR_' TarVar{1}])));
                
                %__________________________________________________________________
            end
            
            %__________________________________________________________________           
            sph = subplot(length(VarList),length(AtlasList),spc);
            hold on; box on; 
            title([VarList2Show{cnt_v} '(' Atlas{1} ')'],'fontsize',fs,'interpreter','latex')
            %__________________________________________________________________
            hndls = ScatterBoxPlots(abs([ndf,hmedf,hcrdf]),'subplot',sph,'color',Col(1:3,:),'pointsize',15,'MedLineColor','k');
            
            sph.XTickLabel=fdMList_labels;
            sph.FontSize = 11;
            sph.XTickLabelRotation=45;
            
            ylabel('Mean FC','fontsize',fs,'interpreter','latex')
            
            ylim(YLims{cnt_v});
            
            Htest(cnt_v,cnt_a) = ttest2(hcrdf,hmedf,'vartype','unequal');
            
            clear ndf hndf Naive Adj
            cnt_a = cnt_a+1;
            spc   = spc+1;
    end

    
    cnt_v = cnt_v+1;
end

set(gcf,'Color','w');

%export_fig(fh,'Global_wCEDen.pdf');