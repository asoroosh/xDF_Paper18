%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
% Soroosh Afyouni, University of Oxford, 2019, 
% srafyouni@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

addpath(genpath('/Users/sorooshafyouni/Home/GitClone/xDF'))

fs = 12;
lw = 1.5;
T_list = [100 200 600 1200]; t_cnt = 1;
SubID = 118528;
%SubID = 135932;

Np = 114;
CG = Np*(Np-1)./2;

nRlz = 2000;

Col=get(groot,'defaultAxesColorOrder');

TVFlag = 'TVOn';

%AlphaS = [0.2 0.5 0.7 1];
AlphaS = ones(1,7);

EstsLables       = {'ME','MEs1','MEt1','MEt2','MEc4','CR','CH','AR1','AR1MC','Fox','Naive'};
EstsLables_names = {'xDF','xDF','xDF','xDF','xDF','BH','Q47','B35','AR1MCPS','G-Q47','Naive'};
OnlyThisMethedos = [5 6:11];
EstsLables       = EstsLables(OnlyThisMethedos);
EstsLables_names = EstsLables_names(OnlyThisMethedos);

for t_cnt = 1:numel(T_list)
    Tp = T_list(t_cnt);
    disp(['T: ' num2str(Tp)])
    for i = 1:nRlz
        if ~mod(i,100); disp(num2str(i)); end; 
            %WhereFrom = ['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FPR/R_' TVFlag '_8-9/Sen_Spec_t' num2str(T_list(t_cnt)) '_' num2str(SubID) '_' num2str(i) '_' TVFlag '_r9.5-8.5_JustZs.mat'];
            %WhereFrom = ['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FPR/R_' TVFlag '_8-9/Sen_Spec_t' num2str(T_list(t_cnt)) '_' num2str(SubID) '_' num2str(i) '_' TVFlag '_r9.5-8.5_JustZs.mat'];
            %WhereFrom = ['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FPR/R_' TVFlag '_2-9/Sen_Spec_t' num2str(T_list(t_cnt)) '_' num2str(SubID) '_' num2str(i) '_' TVFlag '_r9-2_JustZs.mat'];
            WhereFrom = ['/Users/sorooshafyouni/Home/BCF/BCFAnal/Sims/MonsterEquation/FPR/R_' TVFlag '_9-2/Sen_Spec_t' num2str(T_list(t_cnt)) '_' num2str(SubID) '_' num2str(i) '_' TVFlag '_r9-2_JustZs.mat'];
            if ~exist(WhereFrom,'file')
                    disp(['Missing Realisation:' num2str(i)])
                    continue;
            end

            R=load(WhereFrom);
            for m = 1:numel(EstsLables)
                z_raw  = R.Z(:,:,OnlyThisMethedos(m));
                z_raw(1:Np+1:end) = 0;

                p_raw = 2 .* normcdf(-abs(z_raw));

                zt_tmp = fdr_bh(p_raw);

                zGT    = atanh(R.rCovMat).*sqrt(Tp-3);

                POS = numel(find(zGT))./2; %number of all TRUE edges in GROUND TRUTH
                NEG = CG-POS; %number of all Negative edges in the Ground Truth

                TP_tmp = numel(find( zt_tmp .* zGT))./2;
                FN_tmp = numel(find(~zt_tmp .* zGT))./2;

                FP_tmp = numel(find( zt_tmp .* ~zGT))./2;
                TN_tmp = NEG-FP_tmp;

                Sen_tmp = TP_tmp./(TP_tmp+FN_tmp);
                Spc_tmp = TN_tmp./(TN_tmp+FP_tmp);
                Acc_tmp = (TP_tmp+TN_tmp)./(TP_tmp+FP_tmp+FN_tmp+TN_tmp);


                Sen(m,t_cnt,i) = Sen_tmp;
                Spc(m,t_cnt,i) = Spc_tmp;
                Acc(m,t_cnt,i) = Acc_tmp;

            end
    end
end

mSpc = (mean(Spc,3).*100)';
mSen = (mean(Sen,3).*100)';
mACC = (mean(Acc,3).*100)';

sSpc = std(Spc,[],3).*100;
sSen = std(Sen,[],3).*100;
sACC = std(Acc,[],3).*100;

YStone = 110;

fh = figure; 

%-----------------------------
sph0 = subplot(3,1,1);
bh0 = bar(mSpc);
ctr = []; ydt = [];
for k1 = 1:size(mSpc,2)
    ctr_tmp = bsxfun(@plus, bh0(1).XData, [bh0(k1).XOffset]');
    ydt_tmp = bh0(k1).YData;

    ctr(k1,:) = ctr_tmp;
    ydt(k1,:) = ydt_tmp;
    
    bh0(k1).FaceColor = Col(k1,:);
    bh0(k1).FaceAlpha = AlphaS(k1);
        
    for m_cnt = 1:numel(T_list)
        txth = text(ctr_tmp(m_cnt),ydt_tmp(m_cnt)+5,num2str(round(ydt_tmp(m_cnt),2)),'fontsize',10);
        set(txth,'rotation',90)
    end
    
    clear *_tmp
end
ylabel('Specificity','fontsize',fs,'Interpreter','latex')
hold on
%errorbar(ctr', ydt', sSpc,'.r')
ylim([10 150])

%sph0.XTickLabel = EstsLables_names;
sph0.XTickLabel = {'','','',''};
sph0.XTickLabelRotation = 45;
%legend({'T=100','T=200','T=600','T=1200'},'location','southeast')
%-----------------------------
sph1 = subplot(3,1,2);
bh1 = bar(mSen);

ctr = []; ydt = [];
for k1 = 1:size(mSen,2)
    ctr_tmp = bsxfun(@plus, bh1(1).XData, [bh1(k1).XOffset]');
    ydt_tmp = bh1(k1).YData;
    
        bh1(k1).FaceColor = Col(k1,:);
        bh1(k1).FaceAlpha = AlphaS(k1);
        
    for m_cnt = 1:numel(T_list)
        txth = text(ctr_tmp(m_cnt),ydt_tmp(m_cnt)+5,num2str(round(ydt_tmp(m_cnt),2)),'fontsize',10);
        set(txth,'rotation',90)
    end

    ctr(k1,:) = ctr_tmp;
    ydt(k1,:) = ydt_tmp;
    
    clear *_tmp
end
hold on
%errorbar(ctr', ydt', sSen,'.r')
ylabel('Sensitivity','fontsize',fs,'Interpreter','latex')
ylim([0 150])

%sph1.XTickLabel = EstsLables_names;
sph1.XTickLabel = {'','','',''};
sph1.XTickLabelRotation = 45;
%legend({'T=100','T=200','T=600','T=1200'},'location','southeast')
%-----------------------------
sph2 = subplot(3,1,3);
bh2 = bar(mACC);

ctr = []; ydt = [];
for k1 = 1:size(mACC,2) %this is around sample size
    ctr_tmp = bsxfun(@plus, bh2(1).XData, [bh2(k1).XOffset]');
    ydt_tmp = bh2(k1).YData;
    
    %for ii = 1:4; 
        bh2(k1).FaceColor = Col(k1,:);
        bh2(k1).FaceAlpha = AlphaS(k1);
    %end;
        
    for m_cnt = 1:numel(T_list) %this is around methods + Sample Sizes
        txth = text(ctr_tmp(m_cnt),ydt_tmp(m_cnt)+5,num2str(round(ydt_tmp(m_cnt),2)),'fontsize',10);
        set(txth,'rotation',90)
    end
    
    ctr(k1,:) = ctr_tmp;
    ydt(k1,:) = ydt_tmp;
    
    clear *_tmp
end
hold on
%errorbar(ctr', ydt', sACC,'.r')
ylabel('Accuracy','fontsize',fs,'Interpreter','latex')
ylim([10 150])
%sph2.XTickLabel = EstsLables_names;
sph2.XTickLabel = {'','','',''};
sph2.XTickLabelRotation = 45;
set(fh,'color','w')
%legend({'T=100','T=200','T=600','T=1200'},'location','southeast')

legend(EstsLables_names)

%export_fig(fh,'SenSpc.pdf')
