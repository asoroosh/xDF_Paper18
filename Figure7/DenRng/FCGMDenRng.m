%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
% Soroosh Afyouni, University of Oxford, 2019, 
% srafyouni@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Readings%--------------------------------------%--------------------------------------%-----------------------------
load('/home/wmrnaq/HCP/FC/S/HCP_100Unrel_SubList.mat')
SubList = HCP_10Unrel_SubList;

SubID = str2double(getenv('SGE_TASK_ID'));

Sub = SubList{SubID};

%Initialisations%--------------------------------------%--------------------------------------%----------------------

addpath(genpath('~/bin/xDF'))
addpath(genpath('~/bin/2017_01_15_BCT'))


%----------------------%----------------------
DenRng = 0.05:0.05:0.40;
tarden = DenRng(DenID);
disp(['Density is: ' num2str(tarden)])
%----------------------%----------------------

%if strcmp(Atlas,'Gordon')
%	tarden = 0.16;
%elseif strcmp(Atlas,'Power')
%	tarden = 0.15;
%elseif strcmp(Atlas,'Yeo')
%	tarden = 0.23;
%elseif strcmp(Atlas,'ICA200')
%	tarden = 0.1833; %MEs, Naive, CR: [0.18 0.19 0.18]
%elseif strcmp(Atlas,'CC200')
%	tarden = 0.2067; %MEs, Naive, CR: [0.20 0.22 0.2] 
%else
%	error('What?!')
%end

TVFlag = 'TVOff';

from_rt = ['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FCMats/' Atlas '/' PP '/'];
load([from_rt '/Mats_' Atlas '_' PP '_' num2str(Sub) '_' TVFlag '.mat'],'zMats','Info');

save_rt = ['/storage/essicd/data/HCP/Soroosh/HCP_100UR_FCGM/DenRng/' Atlas '/' PP '/' num2str(tarden*100) '/'];

if exist(save_rt,'dir')~=7; mkdir(save_rt); end

MethodList = {'Naive','MEs','CR'}; %MEs corresponds to the ME-Shrink-1

for mthd = MethodList
	if strcmp(mthd{1},'Naive')
		wmat = threshold_proportional(zMats.Naive,tarden);
		bmat = wmat;
		bmat(bmat>0) = 1; 
	elseif strcmp(mthd{1},'MEs')
		wmat = threshold_proportional(zMats.MEs,tarden);
		bmat = wmat;
                bmat(bmat>0) = 1;
        elseif strcmp(mthd{1},'CR')
                wmat = threshold_proportional(zMats.CR,tarden);
                bmat = wmat;
                bmat(bmat>0) = 1;
	end
	%--------------------------------------
%	T = size(bmat,1);	
	N = size(bmat,2)
	
	%--------------------------------------
	disp('- Dgr')
	%DENSITY
	b.Den = density_und(bmat)
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
	save([save_rt '/GM_' Atlas '_' PP '_' mthd{1} '_' num2str(Sub) '_Den' num2str(tarden*100) '_' TVFlag '.mat'],'Info','b','w');

	clear b w
end	
disp('DONE!')
