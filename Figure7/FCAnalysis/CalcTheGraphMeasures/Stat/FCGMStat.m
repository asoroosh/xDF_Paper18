%Readings%--------------------------------------%--------------------------------------%-----------------------------
load('/home/wmrnaq/HCP/FC/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;

SubID = str2double(getenv('SGE_TASK_ID'));

Sub = SubList{SubID};

%Initialisations%--------------------------------------%--------------------------------------%----------------------
addpath(genpath('~/bin/DVARS'))
addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/2017_01_15_BCT'))

TVFlag = 'TVOff';

desrng  = 0.05:0.01:0.50 ;
from_rt = ['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FCMats/' Atlas '/' PP '/'];
load([from_rt '/Mats_' Atlas '_' PP '_' num2str(Sub) '_' TVFlag '.mat'],'twMats','tbMats','Info');

save_rt = ['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FCGM/Stat/' Atlas '/' PP '/'];

if exist(save_rt,'dir')~=7; mkdir(save_rt); end

MethodList = {'Naive','MEs','CR'}; %MEs is for Monster Equation, shrink 1.

for mthd = MethodList
	if strcmp(mthd{1},'Naive')
		bmat = tbMats.Naive_fdr;
		wmat = twMats.Naive_fdr; 
	elseif strcmp(mthd{1},'MEs')
		bmat = tbMats.MEs_fdr;
		wmat = twMats.MEs_fdr;
        elseif strcmp(mthd{1},'CR')
                bmat = tbMats.CR_fdr;
                wmat = twMats.CR_fdr;
	end
	%--------------------------------------
	T = size(bmat,1);	
	N = size(bmat,2)
	%--------------------------------------
	disp('- Dgr')
	%DENSITY
	b.Den = density_und(bmat);
	w.Den = density_und(wmat);

	%DEGREE
	b.Dgr = degrees_und(bmat);
	w.Dgr = strengths_und(wmat);
	
	disp('- Btwn')
	%BTWNSS CENTRALITY
	b.BC = betweenness_bin(bmat);
	w.BC = betweenness_wei(weight_conversion(wmat,'lengths'))./((N-1)*(N-2));

	disp('- Eff - local')
	%LOCAL EFFICIENCY
	b.LE = efficiency_bin(bmat,1);
        w.LE = efficiency_wei(wmat,1);

	disp('- Eff - global')
	%GLOBAL EFFICIENCY
	b.GE = efficiency_bin(bmat);
	w.GE = efficiency_wei(wmat);
	
	%--------------------------------------%--------------------------------------
	save([save_rt '/GM_' Atlas '_' PP '_' mthd{1} '_' num2str(Sub) '_' TVFlag '.mat'],'Info','b','w');

	clear b w
end	
disp('DONE!')
