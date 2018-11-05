clear

fs= 12;

addpath('~/bin/HetBiv')

SubID = 118528;
%SubID = 135932;
TVFlag = 'TVOn';
T_list = [100 200 600 1200];
Np = 114;
CG = Np*(Np-1)./2;
nRlz = 500;

%thrng = 0:0.001:.2;
thrng = 0:1e-4:0.4;

%Titto = figure;

EstsLables = {'ME','MEs1','MEt1','MEt2','MEc4','CR','CH','AR1','AR1MC','Fox','Naive'};

OnlyThisMethedos = [1:11];

%OnlyThisMethedos = [1 6:11];

EstsLables = EstsLables(OnlyThisMethedos);

for t_cnt = 1:numel(T_list)
    Tp = T_list(t_cnt);
    
    disp(['T: ' num2str(t_cnt)])
    
    ii_cnt = 1;
    for ii=1:nRlz
        if ~mod(ii,100); disp(num2str(ii)); end; 
        
        WhereFrom = ['R/Sen_Spec_t' num2str(T_list(t_cnt)) '_' num2str(SubID) '_' num2str(ii) '_' TVFlag '_r9.5-8.5_JustZs.mat'];
            if ~exist(WhereFrom,'file')
                    disp(['Missing Realisation:' num2str(ii)])
                    continue;
            end

            R=load(WhereFrom);
            
%             Acc(ii_cnt,:,t_cnt) = R.Acc;
%             Sen(ii_cnt,:,t_cnt) = R.Sen;
%             Spc(ii_cnt,:,t_cnt) = R.Spc; 
%             FPR(ii_cnt,:,t_cnt) = R.FP./(R.FP+R.TN); 

            for m = 1:numel(EstsLables)
                z_raw  = R.Z(:,:,OnlyThisMethedos(m));
                z_raw(1:Np+1:end) = 0;
                p_raw = 2 .* normcdf(-abs(z_raw));
                
                zGT    = atanh(R.rCovMat).*sqrt(Tp-3);
                th_cnt = 1; 
                for th = thrng

                    %zt_tmp = fdr_bh(p_raw,th);
                    zt_tmp = p_raw<th;
                    zt_tmp(1:Np+1:end) = 0;
                    
                    POS = numel(find(zGT))./2; %number of all TRUE edges in GROUND TRUTH
                    NEG = CG-POS; %number of all Negative edges in the Ground Truth

                    TP_tmp = numel(find( zt_tmp .* zGT))./2;
                    FN_tmp = numel(find(~zt_tmp .* zGT))./2;

                    FP_tmp = numel(find( zt_tmp .* ~zGT))./2;
                    TN_tmp = NEG-FP_tmp;

                    Sen_tmp = TP_tmp./(TP_tmp+FN_tmp);
                    Spc_tmp = TN_tmp./(TN_tmp+FP_tmp);
                    Acc_tmp = (TP_tmp+TN_tmp)./(TP_tmp+FP_tmp+FN_tmp+TN_tmp);
                    
                    Sen_thtmp(th_cnt) = Sen_tmp;
                    Spc_thtmp(th_cnt) = Spc_tmp;

                    Sen(th_cnt,m,t_cnt,ii) = Sen_tmp;
                    Spc(th_cnt,m,t_cnt,ii) = Spc_tmp;
                    Acc(th_cnt,m,t_cnt,ii) = Acc_tmp;

                    clear *_tmp
                    th_cnt = th_cnt+1;
                end
                FPR = 1-Spc_thtmp;
                FPR_LetMeIn = find(FPR<0.1);
                if ~FPR_LetMeIn
                    continue;
                else
                    AUCVal(m,t_cnt,ii) = trapz(FPR(FPR_LetMeIn),Sen_thtmp(FPR_LetMeIn),2);
                end
                clear *_thtmp
            end
            clear R
            ii_cnt = ii_cnt + 1;
    end
end

save(['SenSpcAcc_AUC' TVFlag '_' num2str(nRlz) '.mat'],'Sen','Spc','Acc','AUCVal','thrng')

