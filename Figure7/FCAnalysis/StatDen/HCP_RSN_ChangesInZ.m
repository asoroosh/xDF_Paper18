clear 
%_____________________________________________________
load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FCAnal/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;
SubList(3) = [];

VarList={'Dgr','BC','LE'}; %,'C'};
VarList2Show={'\textbf{Weighted Degree}','\textbf{Betweenness}','\textbf{Local Efficiency}'};
%_____________________________________________________

AtlasList={'Power','Yeo','Gordon'};%,'Gordon'};
GSRList={'GSR'};%,'NoGSR'};
fdMList={'Naive','MEs'};

PP = 'FPP';

thrM_L={'fdr'};

addpath(genpath('/Users/sorooshafyouni/Home/matlab/Ext/DrawingFuncs'))

load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Yeo17_Networks_FuncBlkAss.mat'     ,'YeoZ','YeoU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Power2011_Networks_FuncBlkAss.mat' ,'PwrZ','PwrU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/Gordon_13Networks_FuncBlkAss.mat'  ,'GordonZ','GordonU')
load('/Users/sorooshafyouni/Home/HCP_Scripts/Scripts/HCP_ERMM/FuncNetAss/AAL2_Networks_FuncBlkAss.mat'      ,'AALZ','AALU')

Z={PwrZ;YeoZ;GordonZ};
U={PwrU;YeoU;GordonU};

cnt_v=1; sp_cnt=1; 
for TarVar = VarList
    cnt_a=1;
    for Atlas = AtlasList
        
            rsnidx = Z{cnt_a};
            n_rsn = max(rsnidx);
        
            disp(['A: ' Atlas{1}])
            for s_cnt = 1:numel(SubList)
                GMDir=['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FCAnal/R_GM/Stat/' Atlas{1} '/' PP '/'];

                Naive = load([GMDir 'GM_' Atlas{1} '_' PP '_' fdMList{1} '_' SubList{s_cnt} '.mat']);
                Adj   = load([GMDir 'GM_' Atlas{1} '_' PP '_' fdMList{2} '_' SubList{s_cnt} '.mat']);
                
                ndf  = eval(['Naive.w.' TarVar{1}]);
                hbdf = eval(['Adj.w.'   TarVar{1}]);
                
                for cnt_rsn = 1:n_rsn
                    rsn_naive_gm_tmp = ndf(find(rsnidx  == cnt_rsn));
                    rsn_me_gm_tmp =    hbdf(find(rsnidx == cnt_rsn));
                    
                    rsn_Naive_GM{cnt_v,cnt_a}(cnt_rsn,s_cnt) = mean(rsn_naive_gm_tmp);
                    rsn_ME_GM{cnt_v,cnt_a}(cnt_rsn,s_cnt) = mean(rsn_me_gm_tmp);
                    
                    clear tmp
                end
                %__________________________________________________________________
            end
            %__________________________________________________________________
            
            m_changes_tmp = mean((rsn_ME_GM{cnt_v,cnt_a}-rsn_Naive_GM{cnt_v,cnt_a})./rsn_Naive_GM{cnt_v,cnt_a},2)*100;
            s_changes_tmp = std((rsn_ME_GM{cnt_v,cnt_a}-rsn_Naive_GM{cnt_v,cnt_a})./rsn_Naive_GM{cnt_v,cnt_a},[],2)*100;
            
            spfh = subplot(numel(AtlasList),numel(VarList),sp_cnt);
            hold on; box on; grid on; 
            bph = errorbar(m_changes_tmp,s_changes_tmp,'linestyle','-','linewidth',1.3);
            
            spfh.XTick = 1:max(Z{cnt_a});
            
            xlim([0 max(Z{cnt_a})+1])
            
            if sp_cnt>6
                spfh.XTickLabel = U{cnt_a};
                spfh.XTickLabelRotation = 90;
            end
            
            ChangeRate{cnt_v,cnt_a} = m_changes_tmp;
            
            cnt_a  = cnt_a+1;
            sp_cnt = sp_cnt+1;
    end
    
    cnt_v = cnt_v+1;
end

