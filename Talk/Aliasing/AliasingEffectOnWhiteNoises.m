clear

addpath(genpath('~/Home/GitClone/xDF'))

nRlz = 5000; 

rholist = [0 0.2 0.4 0.6 0.8];
T = 1200; 

for rcnt = 1:numel(rholist)
    rho = rholist(rcnt);
    rho
    for i = 1:nRlz
        ts = mvnrnd([0 0],[1 rho; rho 1],T);
        Rtmp = corr(ts);
        Rtmp = Rtmp(1,2); 
        
        R(rcnt,i)   = Rtmp; 
        Z(rcnt,i)   = atanh(Rtmp);
        sZ(rcnt,i)  = atanh(Rtmp).*sqrt(T-3);
        R2Z(rcnt,i) = Rtmp.*sqrt(T-3);
        
        % Bayley and Hammersely -- Clifford, Richardson & Hemon 
        [Vtmp,ZBHtmp]   = CRBCF(ts',T,[],0);
        VBH(rcnt,i)     = Vtmp(1,2); 
        ZBH(rcnt,i)     = ZBHtmp(1,2); 
        
        [Vtmp,ZBH5tmp]  = CRBCF(ts',T,5,0);
        VBH5(rcnt,i)    = Vtmp(1,2); 
        ZBH5(rcnt,i)    = ZBH5tmp(1,2);
        
        % Bartlett's 1946 --------------------------------------
        [Vtmp,ZB47tmp]  = CheltonBCF(ts',T,[],0); 
        Q47(rcnt,i)     = Vtmp(1,2); 
        ZQ47(rcnt,i)    = ZB47tmp(1,2); 
        
        [Vtmp,ZB475tmp] = CheltonBCF(ts',T,5,0); 
        Q475(rcnt,i)    = Vtmp(1,2);   
        ZQ475(rcnt,i)   = ZB475tmp(1,2);
        
        %[VxDFtmp,ZxDFtmp]= xDF(ts',T,'taper','tukey',sqrt(T));
        [VxDFtmp,ZxDFtmp]= xDF(ts',T,'truncate','adaptive');
        VxDF(rcnt,i) = VxDFtmp(1,2); 
        ZxDF(rcnt,i) = ZxDFtmp.z(1,2);
    end
end


TV = ((1-rholist.^2).^2)./T;

fh0 = figure; 
hold on;  grid on;
plot(var(R')','Marker','o','LineWidth',1.3);
plot(var(Z')','Marker','o','LineWidth',1.3);
fh0.Children.XTick = 1:numel(rholist);
fh0.Children.XTickLabels = {'0','0.2','0.4','0.6','0.8'};

plot(TV,'Color','r','marker','x','linewidth',1.3)

ylabel('$\bf{V}(\rho_{XY})$','Interpreter','latex','fontsize',14)
xlabel('True Correlation Coefficient','fontsize',12)
legend({'Emprical Variance','Stabilised Emprical Variance','Theoritical Variance $N^{-1}(1-\rho_{XY}^2)^2$'},'fontsize',14,'location','northwest','Interpreter','latex')
set(fh0,'Color','w');
ylim([0 1.8e-3])
export_fig(fh0,'VARIANCE_Theoritical&Fisher.pdf')

fh1 = figure; 
hold on;  grid on;
plot(var(R')','Marker','o','LineWidth',1.3);
plot(var(Z')','Marker','o','LineWidth',1.3);
plot(TV,'Color','r','marker','x','linewidth',1.3)
bh = plot([mean(VBH,2) mean(VBH5,2) mean(Q47,2) mean(Q475,2)],'Marker','o','LineWidth',1.3);
fh1.Children.XTick = 1:numel(rholist);
fh1.Children.XTickLabels = {'0','0.2','0.4','0.6','0.8'};

ylabel('$\bf{V}(\rho_{XY})$','Interpreter','latex','fontsize',14)
xlabel('True Correlation Coefficient','fontsize',12)
legend({'Emprical Variance','Stabilised Emprical Variance','Theoritical Variance $N^{-1}(1-\rho_{XY}^2)^2$','BH','BH - Truncated','Q47','Q47 - Fixed Truncation (1/5)'},'fontsize',14,'location','northwest','Interpreter','latex')
set(fh1,'Color','w');
ylim([0 1.8e-3])
export_fig(fh1,'VARIANCE_Theoritical&Fisher&Else.pdf')

fh2 = figure; 
hold on;  grid on;
plot(var(R')','Marker','o','LineWidth',1.3);
plot(var(Z')','Marker','o','LineWidth',1.3);
plot(TV,'Color','r','marker','x','linewidth',1.3)
bh = plot([mean(VBH,2) mean(VBH5,2) mean(Q47,2) mean(Q475,2) mean(VxDF,2)],'Marker','o','LineWidth',1.3);
fh2.Children.XTick = 1:numel(rholist);
fh2.Children.XTickLabels = {'0','0.2','0.4','0.6','0.8'};

ylabel('$\bf{V}(\rho_{XY})$','Interpreter','latex','fontsize',14)
xlabel('True Correlation Coefficient','fontsize',12)
legend({'Emprical Variance','Stabilised Emprical Variance','Theoritical Variance $N^{-1}(1-\rho_{XY}^2)^2$','BH','BH - Truncated','Q47','Q47 - Fixed Truncation (1/5)','xDF - Adaptive Truncation'},'fontsize',14,'location','northwest','Interpreter','latex')
set(fh2,'Color','w');
ylim([0 1.8e-3])
export_fig(fh2,'VARIANCE_Theoritical&Fisher&Else&xDF.pdf')

%%%%%----------------- Z-scores

fh0 = figure; 
hold on;  grid on;
plot(var(R2Z,0,2),'Marker','o','LineWidth',1.3);
plot(var(sZ,0,2),'Marker','o','LineWidth',1.3);
fh0.Children.XTick = 1:numel(rholist);
fh0.Children.XTickLabels = {'0','0.2','0.4','0.6','0.8'};

ylabel('$\bf{V}(Z_{XY})$','Interpreter','latex','fontsize',14)
xlabel('True Correlation Coefficient','fontsize',12)
legend({'Empirical Z-scores','Stabilised Emprical Z-scores'},'fontsize',14,'location','southwest','Interpreter','latex')
set(fh0,'Color','w');
%ylim([0 0.6e-3])
export_fig(fh0,'Z_Theoritical&Fisher.pdf')

fh1 = figure; 
hold on;  grid on;
plot(var(R2Z')','Marker','o','LineWidth',1.3);
plot(var(sZ')','Marker','o','LineWidth',1.3);
bh = plot([var(ZBH,0,2) var(ZBH5,0,2) var(ZQ47,0,2) var(ZQ475,0,2)],'Marker','o','LineWidth',1.3);
fh1.Children.XTick = 1:numel(rholist);
fh1.Children.XTickLabels = {'0','0.2','0.4','0.6','0.8'};

ylabel('$\bf{V}(Z_{XY})$','Interpreter','latex','fontsize',14)
xlabel('True Correlation Coefficient','fontsize',12)
legend({'Emprical Z-scores','Stabilised Emprical Z-scores','BH','BH - Fixed Truncation (1/5)','Q47','Q47 - Fixed Truncation (1/5)'},'fontsize',14,'location','southwest','Interpreter','latex')
set(fh1,'Color','w');
%ylim([0 0.6e-3])
export_fig(fh1,'Z_Theoritical&Fisher&Else.pdf')

 fh2 = figure; 
 hold on;  grid on;
 plot(var(R2Z')','Marker','o','LineWidth',1.3);
 plot(var(sZ')','Marker','o','LineWidth',1.3);
 bh = plot([var(ZBH,0,2) var(ZBH5,0,2) var(ZQ47,0,2) var(ZQ475,0,2) var(ZxDF,0,2)],'Marker','o','LineWidth',1.3);
 fh2.Children.XTick = 1:numel(rholist);
 fh2.Children.XTickLabels = {'0','0.2','0.4','0.6','0.8'};
 
 ylabel('$\bf{V}(Z_{XY})$','Interpreter','latex','fontsize',14)
 xlabel('True Correlation Coefficient','fontsize',12)
 legend({'Emprical Z-scores','Stabilised Emprical Z-scores','BH','BH - Fixed Truncation (1/5)','Q47','Q47 - Fixed Truncation (1/5)','xDF - Adaptive Truncation'},'fontsize',14,'location','southwest','Interpreter','latex')
 set(fh2,'Color','w');
 %ylim([0 0.6e-3])
 export_fig(fh2,'Z_Theoritical&Fisher&Else&xDF.pdf')