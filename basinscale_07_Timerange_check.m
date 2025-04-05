%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to verify the optimal layer view of striations by using idea isopycnal
% model "Aronnax".
% check the flow field after highpass filtering
%  ------------------------------------------------------------------------
% by zhangyu 2020526
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  initializaiotn
clear all
addpath(genpath('/home/dell/00_MatlabProgram/'));
addpath(genpath('/home/dell/Striations/'));

pth0 = '/home/dell/Striations/striations_3layers_basinscale/';

R = 6371*1.e3;      % radius of the earth, unit: m

% basinscale_3layer_s0
xlen0 = [2880, 3840, 5760, 7680];
ylen0 = [3840, 3840, 3840, 3840];
nx0 =   [ 97, 129, 193, 257];
ny0 =   [129, 129, 129, 129];

% basinscale_3layer_re15km
xlen1 = [ 960, 1920, 2880, 3840, 4800];
ylen1 = [3840, 3840, 3840, 3840, 3840];
nx1 = [ 64, 128, 192, 256, 320];
ny1 = [256, 256, 256, 256, 256];

yt = 1:10:150;

%% read and plot >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for s = 5 % loop of scenarios

xlen = xlen1(s)*1.e3;
ylen = ylen1(s)*1.e3;
nx = nx1(s);
ny = ny1(s);
dx = xlen/nx;
dy = ylen/ny;
re = dx/1.e3/100;

dpth = [pth0,'s1', num2str(s),'/'];

figure('color','w','position',[100 100 1500 800])
aa = [12,13,14,15];
for t = 1:length(aa)  % loop of timerange
	y1 = yt(aa(t)) 
	y2 = yt(aa(t)) +9
load([dpth,'output_eta_mean_',num2str(y1),'-',num2str(y2),'.mat']);
load([dpth,'output_h_mean_',num2str(y1),'-',num2str(y2),'.mat']);
load([dpth,'output_u_mean_',num2str(y1),'-',num2str(y2),'.mat']);
load([dpth,'output_v_mean_',num2str(y1),'-',num2str(y2),'.mat']);


h(1,:,:) = h(1,:,:) + eta;

% u_0 and u_str ================================================
%figure('color','w','position',[100 100 1500 800])

for i = 1:3    % loop of layers
    tmp_u = squeeze(u(i,2:end,:));size(tmp_u)
    tmp_v = squeeze(v(i,:,2:end));size(tmp_v)
    [X,Y] = meshgrid(dx:dx:xlen,dy:dy:ylen);
    
    % streamfunction
    [PHI,PSI] = flowfun(tmp_u,tmp_v);

    %{
    % u_0 ------------------------------------
    subplot(2,3,i)
    hold on
    contourf(X,Y,tmp_u',100,'linestyle','none');
    colorbar;
    %cmocean('balance','pivot',0)
    colormap(create_colors(100-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]));
    caxis([-.1 .1]);    
    
    % steamlines
    % contour(X,Y,PSI,10,'-k','linewidth',1.2);
    
    title(['S',num2str(s),': u_0 on layer ',num2str(i)],'fontweight','bold');
    xlabel('x (1e3 km)');ylabel('y (1e3 km)');
    set(gca,'xtick',1.e6:1.e6:8.e6,'xticklabel',{'1','2','3','4','5','6','7','8'}...
        ,'ytick',1.e6:1.e6:8.e6,'yticklabel',{'1','2','3','4','5','6','7','8'});
    set(gca,'fontsize',12,'box','off','tickdir','out','linewidth',1);
    %}

    % u_str -------------------------------------
    cutoff_length = 400;
    tmp_us = striations_reveal(tmp_u,re,cutoff_length);
    
    subplot(3, 4, t+(i-1)*4)
    vel_contours = [-100,-30:10:-10,-9:-2,-1.5,-1:0.2:1,1.5,2:9,10:10:30,100]/10;
    contourf(X,Y,tmp_us'*100,vel_contours,'linestyle','none');
    
    if t==4 & i==2
    cb = colorbar('position',[0.92 0.15 0.02 0.7]);
    set(get(cb,'title'),'string','cm s^{-1}');
    set(cb,'tickdir','out')
    end
    
    colormap(create_colors(100-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]));
    caxis([-.5 .5]);
    
    title([num2str(y1),'-',num2str(y2),'th : u_{s} on layer ',num2str(i)],'fontweight','bold');
    
    if i == 3;xlabel('x (km)');end
    if s == 1;ylabel('y (km)');end
    
    if s ==1
            set(gca,'xtick',2.e5:2.e5:1.e6,'xticklabel',{'200','400','600','800','1000'});
    elseif s ==2
	    set(gca,'xtick',5.e5:5.e5:2.e6,'xticklabel',{'500','1000','1500','2000'});
    else
	    set(gca,'xtick',1.e6:1.e6:5.e6,'xticklabel',{'1000','2000','3000','4000','5000'});
    end
    set(gca,'ytick',1.e6:1.e6:4.e6,'yticklabel',{'1000','2000','3000','4000'});
    set(gca,'fontsize',12,'box','off','tickdir','out','linewidth',1);

end % end of for i = 

end % end of for t = 

end % end of for s =
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
