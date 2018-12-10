clear

load('/home/wmrnaq/HCP/FC/S/HCP_100Unrel_SubList.mat')
SubList=HCP_10Unrel_SubList;

addpath('/home/wmrnaq/bin/HetBiv/StatThresholding')
addpath('/home/wmrnaq/bin/2017_01_15_BCT')

SubList(3)=[];

Atlas = 'CC200';

densrng = 0.01:0.02:0.50;

for i=1:99
    
    disp(['Sub: ' num2str(i) ' - ' SubList{i}])
    
    clear mat
    load(['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FCMats/' Atlas '/FPP/Mats_' Atlas '_FPP_' SubList{i} '_TVOn.mat'],'zMats')
    mat = zMats.CR;
    
    [REST_CE(i),REST_CEDns(i)] = CostEff_bin(mat,densrng);    
end

save(['HCP_100UR_' Atlas '_CE.mat'],'REST_CE','REST_CEDns')
    
    
