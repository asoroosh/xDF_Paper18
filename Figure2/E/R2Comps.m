
clear

lw =1.3;
fs = 12;

AtlasList = {'Gordon','ICA200','Power','Yeo','MMP'};

Col = get(groot,'defaultAxesColorOrder');

a_cnt = 1; 
for Atlas = AtlasList

    load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/ROIsVSCorrLeng/R/HCP100UR/' Atlas{1} '_GSR_FPP_CorrLeng_ROIWise.mat'])
    
    [N,nSub] = size(ACL);
    
    Y = ACL'; 
    
    SStot = sum(sum((Y - mean(Y(:))).^2));
    SSerr = sum(sum((Y - repmat(mean(Y),[nSub 1])).^2));
    InterSubR2(a_cnt) = (SStot-SSerr)/SStot;

    
    Y = Y'; 
    
    SStot = sum(sum((Y - mean(Y(:))).^2));
    SSerr = sum(sum((Y - repmat(mean(Y),[N 1])).^2));
    InterNodeR2(a_cnt) = (SStot-SSerr)/SStot;
    
    a_cnt = a_cnt+1; 
end

fh = figure('position',[50,500,300,250]); 
hold on; box on; grid on; 
bh = bar(([InterSubR2;InterNodeR2]').*100);
ylabel('Variance explained in ACL profile (R$^2$)','FontSize',fs,'Interpreter','latex')
xlabel('Parcellation Scheme','FontSize',fs,'Interpreter','latex')
fh.Children.XTick = 1:numel(AtlasList);
fh.Children.XTickLabel = AtlasList;
fh.Children.XTickLabelRotation = 45; 
bh(1).FaceColor = Col(1,:);
bh(2).FaceColor = Col(2,:);

legend({'Inter-Subject','Inter-Node'},'FontSize',fs)
set(fh,'Color','w');
export_fig(fh,'Figs/ACprofiles.pdf')
