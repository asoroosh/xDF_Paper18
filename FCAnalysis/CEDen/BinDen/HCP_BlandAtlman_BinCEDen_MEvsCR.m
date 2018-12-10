clear 
%_____________________________________________________
load('/Users/sorooshafyouni/Home/BCF/BCFAnal/FCAnal/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;
SubList(3) = [];

VarList={'Dgr','BC','LE'}; %,'C'};
VarList2Show={'\textbf{Weighted Degree}','\textbf{Betweenness}','\textbf{Local Efficiency}'};
%_____________________________________________________

AtlasList={'Power','Yeo','Gordon'};%,'Gordon'};
GSRList={'GSR'};%,'NoGSR'};
fdMList={'MEs','CR'};

PP = 'FPP';

thrM_L={'fdr'};

addpath(genpath('/Users/sorooshafyouni/Home/matlab/Ext/DrawingFuncs'))

load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Yeo17_Networks_FuncBlkAss.mat'     ,'YeoZ','YeoU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Power2011_Networks_FuncBlkAss.mat' ,'PwrZ','PwrU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Gordon_13Networks_FuncBlkAss.mat'  ,'GordonZ','GordonU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/AAL2_Networks_FuncBlkAss.mat'      ,'AALZ','AALU')

Z={PwrZ;YeoZ;GordonZ};
U={PwrU;YeoU;GordonU};

cnt_g=1; spc=1;
    fh = figure('position',[50,500,620,500]);
    hold on; box on;    
cnt_v=1;
for TarVar = VarList
    cnt_a=1;
    for Atlas = AtlasList
            disp(['A: ' Atlas{1}])
            for s_cnt = 1:numel(SubList)
                GMDir=['/Users/sorooshafyouni/Home/BCF/BCFAnal/FCAnal/R_GM/CE/' Atlas{1} '/' PP '/'];

                MEs = load([GMDir 'GM_' Atlas{1} '_' PP '_' fdMList{1} '_' SubList{s_cnt} '_TVOff.mat']);
                CR   = load([GMDir 'GM_' Atlas{1} '_' PP '_' fdMList{2} '_' SubList{s_cnt} '_TVOff.mat']);
                %__________________________________________________________________           
                spfh = subplot(length(VarList),length(AtlasList),spc);
                hold on; box on; 
                title([VarList2Show{cnt_v} '(' Atlas{1} ')'],'fontsize',12,'interpreter','latex')
                %__________________________________________________________________
                ME_df(:,s_cnt)  = eval(['MEs.b.' TarVar{1}]);
                CR_df(:,s_cnt)  = eval(['CR.b.'   TarVar{1}]);
                %__________________________________________________________________
            end
            [~,p]=ttest2(ME_df',CR_df','vartype','unequal'); %Welch Test
            h=fdr_bh(p);
            %h = p<0.05;
            DiffRate=sum(h)./numel(h);
            %__________________________________________________________________
            size(ME_df)
            [HDI]=HCP_BlandAltman_StatDen(ME_df,CR_df,'assignments',Z{cnt_a},'assignmentslabels',U{cnt_a},'test',h);
            %__________________________________________________________________
            xlim0=xlim; ylim0=ylim;
            text((xlim0(1)+(xlim0(2)-xlim0(1))*1/6),(ylim0(1)+(ylim0(2)-ylim0(1))*0.7/6),['B=' num2str(round(HDI,2))]);
            legend(gca,'off')
            %--
            XTickLabel = get(spfh,'XTick');
            set(spfh,'XTickLabel',num2str(XTickLabel'))
            YTickLabel = get(spfh,'YTick');
            set(spfh,'YTickLabel',num2str(YTickLabel'))             
            %__________________________________________________________________
            xlabel(' ')
            ylabel(' ')
            %if ismember(spc,[2,5])
                xlabel('$<$BH$,$xDF$>$','fontsize',13,'interpreter','latex')
            %end
            %if ismember(spc,[1,4])
                ylabel('xDF$-$BH','fontsize',13,'interpreter','latex')
            %end
            
            AllDiffRate(cnt_v,cnt_a) = DiffRate;
            
            cnt_a = cnt_a+1;
            spc   = spc+1;
            clear ME_* CR_*
    end
    cnt_v = cnt_v+1;
end

AllDiffRate*100

set(gcf,'color','w');

writetable(array2table(AllDiffRate*100,'VariableNames',AtlasList,'RowNames',VarList),'Figs/xDFvsBH.csv')

export_fig(fh,'Figs/Local_BinCEDen_MEvsCR.pdf')

