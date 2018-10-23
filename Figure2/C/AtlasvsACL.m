clear

SiteList={'HCP'};
AtlasList={'ICA20050','Gordon','Yeo','MMP'};
GSRList={'GSR'};

load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/ROIsVSCorrLeng/R/ROISizes_Yeo17_AAL2_AAL_CC200_Power2011_ICA2002550.mat']);

Col       = get(groot,'defaultAxesColorOrder');
HCP_AtlasCol  = [Col(5,:);Col(4,:);Col(1,:);Col(3,:)]; %Col(6,:) & Col(7,:) for PCP


sp1=2;
sp2=numel(AtlasList)./2;

for Site=SiteList
    gsr_cnt = 1;
    for GSR=GSRList

        corrfigs = figure('position',[50,500,820,400]);
        hold on; box on;  

%         sizefig = figure('position',[50,500,800,400]);
        
        a_cnt=1;
        for Atlas=AtlasList
            if strcmp(Atlas{1},'MMP')
                load('/Users/sorooshafyouni/Home/BCF/BCFAnal/ROIsVSCorrLeng/R/MMP_vrtxNumber.mat','ROISize')
                ROI_sz=ROISize;
            else
                ROI_sz=eval([Atlas{1} '_nvox_r']);
            end
                if sum(ismember(Atlas{1},'ICA100'))>4
                    Atlas{1} = 'ICA200';
                end
                
                if strcmp(Atlas{1},'MMP')
                    load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/ROIsVSCorrLeng/R/HCP100UR/' Atlas{1} '_NoGSR_FPP_CorrLeng_ROIWise.mat']);
                else
                    load(['/Users/sorooshafyouni/Home/BCF/BCFAnal/ROIsVSCorrLeng/R/HCP100UR/' Atlas{1} '_' GSR{1} '_FPP_CorrLeng_ROIWise.mat']);
                end
                
                nn=size(ACL,1);
                ns=size(ACL,2);
                
                figure(corrfigs)
                subplot(sp1,sp2,a_cnt)
                title([Site{1} ' ' Atlas{1}],'interpreter','latex','fontsize',13)
                hold on; box on; grid on;
                if strcmp(Atlas{1},'MMP')
                    scatter(ROI_sz.^(1/2),mean(ACL,2),'MarkerFaceColor',HCP_AtlasCol(a_cnt,:),'MarkerEdgeColor',HCP_AtlasCol(a_cnt,:));
                    ylabel('Auto-correlation Length ($\tau$)','interpreter','latex','fontsize',12)
                    xlabel('$\sqrt[2]{ROI sizes}~~(vertex)$','interpreter','latex','fontsize',12)
                else
                    scatter(ROI_sz.^(1/3),mean(ACL,2),'MarkerFaceColor',HCP_AtlasCol(a_cnt,:),'MarkerEdgeColor',HCP_AtlasCol(a_cnt,:));
                    ylabel('Auto-correlation Length ($\tau$)','interpreter','latex','fontsize',12)
                    xlabel('$\sqrt[3]{ROI sizes}~~(voxel)$','interpreter','latex','fontsize',12)
                end
                lsline
                %plot(0:20,0:20,'r-.')
       
                ylim([0 8])
                if strcmp(Atlas{1},'MMP')
                    xlim([0 35])
                else
                    xlim([0 15])
                end
                  
                
                [r,p]=corr(ROI_sz',mean(ACL,2));
                
                if ~round(p,3)
                    p = 0.001;
                end
                
                xlim0=xlim; ylim0=ylim;
                text((xlim0(1)+(xlim0(2)-xlim0(1))*4/6),(ylim0(1)+(ylim0(2)-ylim0(1))*2.5/6),  ['r=' num2str(round(r,2))],'FontSize',12);%,'FontWeight','bold');
                text((xlim0(1)+(xlim0(2)-xlim0(1))*4/6),(ylim0(1)+(ylim0(2)-ylim0(1))*2/6),['p<' num2str(round(p,4))],'FontSize',12)%,'FontWeight','bold');
                
                AllACL{a_cnt,gsr_cnt} = ACL;

                ROIs{a_cnt} = ROI_sz;
                
                %clear CorrLeng ROI_sz
                a_cnt=a_cnt+1;
        end
        
        set(gcf,'Color','w')
        
        export_fig(['Figs/CorrBtwn_ROISize_CorrLeng_' Site{1} '_' GSR{1} '.pdf'])
        
        gsr_cnt = gsr_cnt+1;
    end
end
