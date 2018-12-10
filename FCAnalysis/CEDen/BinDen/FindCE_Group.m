clear

load('/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FCAnal/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;
SubList(3) = [];

AtlasList={'Power','Yeo','Gordon'};

PP = 'FPP';

a_cnt = 1;
for Atlas = AtlasList
    for i=1:numel(SubList)
        WhereFrom = ['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FCAnal/R_mats/' Atlas{1} '/' PP '/'];
        load([WhereFrom 'Mats_' Atlas{1} '_' PP '_' num2str(SubList{i}) '.mat'],'twMats','tbMats')

        disp(['Subject: ' num2str(i)])

        % DENSITY ------------------------------------------------------------
        CE_Den(i,1,a_cnt) = density_und(twMats.MEs_CE);
        CE_Den(i,2,a_cnt) = density_und(twMats.CR_CE);
        CE_Den(i,3,a_cnt) = density_und(twMats.Naive_CE);
    end
    a_cnt = a_cnt+1; 
end