clear

fs = 12; 

load('R/HCP_InterSub_DyConn_YEO.mat')

%[~,~,kstmp]=kstest(MthdNow);

Col = get(groot,'defaultAxesColorOrder');

MethodList = {'naive','xdff','xdf'};

fh = figure;
subplot(1,2,1)
hold on; box on; grid on;
title('')
for i = 1:numel(MethodList)
    tmp  = eval(['mat_' MethodList{i}]);
    QQ = qqplot(tmp(1:1000:end)); 
    QQ(1).Marker = '.';
    QQ(1).MarkerEdgeColor = Col(i,:);
    
    QQ(2).Visible = 'off';
    QQ(3).Visible = 'off';
    
    [~,~,KS(i)] = kstest(tmp);
end

rfline = refline(1,0); 
rfline.LineWidth = 1.3;
rfline.Color = 'k';

set(fh,'Color','w');

sph = subplot(1,2,2);
hold on; box on; grid off; 
for i=1:numel(MethodList)
    tq_tmp = -log10(KS(i));
    barfh = bar(i,tq_tmp,'FaceColor',Col(i,:),'FaceAlpha',0.5);
    txth = text(i,tq_tmp+0.01,num2str(round(tq_tmp,3)));
    set(txth,'rotation',90)
end
ylabel(['-log10(KS)'],'fontsize',fs)

ylim([0 2.1])
title('')
sph.XTick = 1:numel(MethodList);
sph.XTickLabel = {'Naive','Window-xDF','xDF'};
sph.XTickLabelRotation = 45;

set(fh,'Color','w');
export_fig(fh,'drsFC_xDF.pdf')