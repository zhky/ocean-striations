%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to verify the optimal layer view of striations by using idea isopycnal
% model "Aronnax".
% check the flow field after highpass filtering
%  ------------------------------------------------------------------------
% by zhangyu 2020526
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
clear all
addpath(genpath('/home/dell/00_MatlabProgram/'));
addpath(genpath('/home/dell/00_MatlabProgram_Striations/'));

re = 2000/128/100;
R = 6371*1.e3;      % radius of the earth, unit: m

xlen1 = [2000*1e3, 4000*1e3, 6000*1e3, 2880*1e3, 3840*1e3, 5760*1e3, 7680*1e3];
ylen1 = [4000*1e3, 4000*1e3, 4000*1e3, 3840*1e3, 3840*1e3, 3840*1e3, 3840*1e3];
nx1 = [128, 256, 384, 96, 129, 193, 257];
ny1 = [256, 256, 256, 129, 129, 129, 129];

%% four scenaros 
for s = 4:7
xlen = xlen1(s);
ylen = ylen1(s);
nx = nx1(s);
ny = ny1(s);
dx = xlen/nx;
dy = ylen/ny;

dpth = ['/home/dell/00_MatlabProgram_Striations/striations_3layers_basinscale/s0',...
    num2str(s),'/'];

load([dpth,'output_eta_mean.mat']);
load([dpth,'output_h_mean.mat']);
load([dpth,'output_u_mean.mat']);
load([dpth,'output_v_mean.mat']);

h(1,:,:) = h(1,:,:) + eta;

% u_0 and u_str ================================================
figure('color','w','position',[100 100 1500 800])

for i = 1:3
    tmp_u = squeeze(u(i,2:end,:));size(tmp_u)
    tmp_v = squeeze(v(i,:,2:end));size(tmp_v)
    [X,Y] = meshgrid(dx:dx:xlen,dy:dy:ylen);
    
    % streamfunction
    [PHI,PSI] = flowfun(tmp_u,tmp_v);

    % u_0
    subplot(2,3,i)
    hold on
    contourf(X,Y,tmp_u',100,'linestyle','none');
    colorbar;
    %cmocean('balance','pivot',0)
    colormap(create_colors(100-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]));
    caxis([-.1 .1]);    
%     switch s
%         case 1; caxis([-1e-8 1e-8]);
%         case 2; caxis([-1e-7 1e-7]);
%         case 3; caxis([-1e-5 1e-5]);
%         case 4; caxis([-1e-5 1e-5]);
%     end
    
%     % steamlines
%     contour(X,Y,PSI,10,'-k','linewidth',1.2);
    
    title(['S',num2str(s),': u_0 on layer ',num2str(i)],'fontweight','bold');
    xlabel('x (1e3 km)');ylabel('y (1e3 km)');
    set(gca,'xtick',1.e6:1.e6:8.e6,'xticklabel',{'1','2','3','4','5','6','7','8'}...
        ,'ytick',1.e6:1.e6:8.e6,'yticklabel',{'1','2','3','4','5','6','7','8'});
    set(gca,'fontsize',12,'box','off','tickdir','out','linewidth',1);

    % u_str
    cutoff_length = 400;
    tmp_us = striations_reveal(tmp_u,re,cutoff_length);
    
    subplot(2,3,i+3)
    contourf(X,Y,tmp_us',100,'linestyle','none');
    colorbar;
    %cmocean('balance','pivot',0);
    
    colormap(create_colors(100-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]));
    caxis([-0.05 0.05]);
%     switch s
%         case 1; caxis([-1e-7 1e-7]);
%         case 2; caxis([-1e-7 1e-7]);
%         case 3; caxis([-1e-5 1e-5]);
%         case 4; caxis([-1e-5 1e-5]);
%     end
   
    title(['S',num2str(s),': u_{s} on layer ',num2str(i)],'fontweight','bold');
    xlabel('x (1e3 km)');ylabel('y (1e3 km)');
    set(gca,'xtick',1.e6:1.e6:8.e6,'xticklabel',{'1','2','3','4','5','6','7','8'}...
        ,'ytick',1.e6:1.e6:8.e6,'yticklabel',{'1','2','3','4','5','6','7','8'});
    set(gca,'fontsize',12,'box','off','tickdir','out','linewidth',1);
end

end

%% streamline =========================================================
% for s = 1:4
%     
%     eval(['load output_0',num2str(s),'_eta_mean.mat']);
%     eval(['load output_0',num2str(s),'_h_mean.mat']);
%     eval(['load output_0',num2str(s),'_u_mean.mat']);
%     eval(['load output_0',num2str(s),'_v_mean.mat']);
% 
%     figure('color','w')
%     for i = 1:3
%         tmp_u = squeeze(u(i,:,2:end));
%         tmp_v = squeeze(v(i,2:end,:));
% 
%         % streamlines
%         subplot(1,3,i)
%         [X,Y] = meshgrid(dx:dx:xlen,dy:dy:ylen);
%         st = streamslice(X,Y,tmp_u,tmp_v);
%         set(st,'linewidth',1.2,'color','b');
%         title(['Streamlines on layer ',num2str(i)]);
%         xlabel('x (km)');ylabel('y (km)');
%         set(gca,'xtick',1.e6:1.e6:5.e6,'xticklabel',{'1000','2000','3000','4000','5000'}...
%             ,'ytick',1.e6:1.e6:5.e6,'yticklabel',{'1000','2000','3000','4000','5000'});
%         set(gca,'fontsize',12,'box','off','tickdir','out','linewidth',1);
%         box on; axis tight;
%     end
% 
% end
