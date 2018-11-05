addpath(genpath('~/bin/DVARS'))
addpath(genpath('~/bin/HetBiv'))
addpath(genpath('~/bin/2017_01_15_BCT'))

n = str2double(getenv('SGE_TASK_ID')); 

rng(n)

save_rt = ['/storage/essicd/data/HCP/Soroosh/DyConn/InterSubMats/' Atlas '/'];
if exist(save_rt,'dir')~=7; mkdir(save_rt); end

load(['/home/wmrnaq/ACAnal/InterSub/InterSubMats/S/MTS_HCP100UR_' Atlas '_GSRed.mat'])

[Nsub,T,N] = size(MTS);

Ntar = 20;

AllPerms = Ntar^2;

WinLngt = 100;
Wratio = T./WinLngt;

EndPoint = T - WinLngt;

TapWind = ones(1,WinLngt);

cnt = 1;

disp(['Intersub scrambling...'])
    w = randi(EndPoint);
    %if ~mod(n,100); disp(['Rlz: ' num2str(n)]); end;
    for i = 1:Ntar
        irandnode = randi(N);
        femalesub = randi(Nsub);
        femalets = MTS(femalesub,:,irandnode)';
        wfemalets = femalets;
        wfemalets = wfemalets(w+(0:WinLngt-1)).*TapWind';
        for j = 1:Ntar
            %Apply the window with multiplication
            jrandnode = randi(N);
            subrng = 1:Nsub;
            subrng(femalesub) = [];
            malesub = datasample(subrng,1,'Replace',false);
            malets  = MTS(malesub,:,jrandnode)';
            wmalets = malets;
            wmalets = wmalets(w+(0:WinLngt-1)).*TapWind';
            
            Mat.mat_naive(i,j) = atanh(corr(wfemalets,wmalets)).*sqrt(WinLngt-3);
                        
            %xDF on the whole time series
            ASAt = xDF([malets,femalets]',T,'taper','shrink');

	    %xDF on the whole time series
            [~,Stat] = xDF([wmalets,wfemalets]',WinLngt,'taper','shrink');
	    Mat.mat_xdff(i,j) = Stat.z.rzf(1,2);            

            %xDF divided by number of windows
            rmat_tmp = corr(wfemalets,wmalets);
            Mat.mat_xdf(i,j) = atanh(rmat_tmp)./(sqrt((ASAt(1,2).*Wratio)./((1-rmat_tmp.^2).^2)));

	    %Prewhitening
	    wts=PreWhitenMe([malets,femalets]',T,'taper','tukey',sqrt(T),'DM','cholesky')';

	    wwts = wts(w+(0:WinLngt-1),:).*repmat(TapWind',[1 2]);
	    mat_tmp=corr(wwts);
	    Mat.mat_pw(i,j) = atanh(mat_tmp(1,2)).*sqrt(WinLngt-3);            

            %RandInd(:,n) = [w irandnode jrandnode malesub femalesub];
            
            disp(num2str(round(cnt/AllPerms,2)))
   
	    cnt = cnt + 1; 
 
        end
    end
    
    Mat.mat_wn = atanh(corr(randn(WinLngt,N))).*sqrt(WinLngt-3);
    Mat.mat_wn(1:N+1:end) = 0;

save([save_rt '/InterSubMats_' Atlas '_' num2str(n) '.mat'],'Mat');

disp('Done!')

