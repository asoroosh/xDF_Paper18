clear

Atlas={'Gordon','Yeo17','Power2011','ICA100','MMP'};
GSR={'GSR'};
%addpath /Users/sorooshafyouni/Home/HCP_Scripts/Scripts/Netmats

addpath(genpath('~/bin/xDF'))
addpath(genpath('~/bin/DVARS'))

load('/home/wmrnaq/HCP/FC/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;

T = 1200;

for a=Atlas
    for pp=PP
            disp([a{1} '_' g{1} ])
            WhereFrom = ['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FC/' a{1} '/' pp{1} '/'];
        for s=1:ns
            load([WhereFrom 'HCP_' pp{1} '_' num2str(SubList{s}) '_ROIsTS.mat'],'mts')
            disp(['Sub: ' num2str(s)])
            
            mts = GSRme(mts,T);
            AC  = AC_fft(mts,T);
            ACL = sum(AC.^2,2);
        end
        clear TS
        save(['R/' a{1} '_' GSR{1} '_' pp{1} '_CorrLeng_ROIWise.mat'],'ACL')
        clear CorrLeng
    end
end