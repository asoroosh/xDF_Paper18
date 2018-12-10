%Readings%--------------------------------------%--------------------------------------%-----------------------------
load('/home/wmrnaq/HCP/FC/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;

SubID = str2double(getenv('SGE_TASK_ID'));

Sub = SubList{SubID};

TVFlag = 'TVOff';

%Initialisations%--------------------------------------%--------------------------------------%----------------------

addpath(genpath('~/bin/DVARS'))
addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/2017_01_15_BCT'))

from_rt = ['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FCMats/' Atlas '/' PP '/'];
load([from_rt '/Mats_' Atlas '_' PP '_' num2str(Sub) '_' TVFlag '.mat'],'zMats','Info');

save_rt = ['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FCGM/Unthr/' Atlas '/' PP '/'];

if exist(save_rt,'dir')~=7; mkdir(save_rt); end

MethodList = {'Naive','MEs','CR'}; %MEs corresponds to the ME-Shrink-1

% I removed the bin versions, as there is no point of dealing with them. 
% Also, I have removed the negatives, because I no one knows how to deal with them in graphg measures. 

for mthd = MethodList
	if strcmp(mthd{1},'Naive')
		wmat = zMats.Naive;
	elseif strcmp(mthd{1},'MEs')
		wmat = zMats.MEs;
        elseif strcmp(mthd{1},'CR')
                wmat = zMats.CR;
	end
	%--------------------------------------
%	T = size(bmat,1);	
	N = size(wmat,2)
	%--------------------------------------
	disp('- Dgr')
	%DENSITY
	w.Den = density_und(wmat);

	%DEGREE
	w.Dgr = strengths_und(wmat);
	
	disp('- Btwn')
	%BTWNSS CENTRALITY
	w.BC = betweenness_wei(weight_conversion(wmat,'lengths'))./((N-1)*(N-2));

	disp('- Eff - local')
	%LOCAL EFFICIENCY
        w.LE = efficiency_wei(wmat,1);

	disp('- Eff - global')
	%GLOBAL EFFICIENCY
	w.GE = efficiency_wei(wmat);
	
	%--------------------------------------%--------------------------------------
	save([save_rt '/GM_' Atlas '_' PP '_' mthd{1} '_' num2str(Sub) '_' TVFlag '.mat'],'Info','w');

	clear b w
end	
disp('DONE!')
