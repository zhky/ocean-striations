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
addpath(genpath('/home/dell/Striations/'));

R = 6371*1.e3;      % radius of the earth, unit: m
cl = [200,300,400];

pth0 = ['/home/dell/Striations/striations_3layers_basinscale/'];

% basinscale_3layer_s0
xlen0 = [2880, 3840, 5760, 7680];
ylen0 = [3840, 3840, 3840, 3840];
nx0 =   [ 97, 129, 193, 257];
ny0 =   [129, 129, 129, 129];

% basinscale_3layer_re15km
xlen1 = [ 960, 1920, 2880, 3840, 4800, 5760];
ylen1 = [3840, 3840, 3840, 3840, 3840, 3840];
nx1 = [ 64, 128, 192, 256, 320, 384];
ny1 = [256, 256, 256, 256, 256, 256];


%% scenaros 
figure('color','w','position',[100 100 1400 800])
for s = 1:6
    s
xlen = xlen1(s)*1.e3;
ylen = ylen1(s)*1.e3;
nx = nx1(s);
ny = ny1(s);
dx = xlen/nx;
dy = ylen/ny;
re = dx/1.e3/100;

xx = dx:dx:xlen;
yy = dy:dy:ylen;

dpth = [pth0,'s1',num2str(s),'/'];

load([dpth,'output_eta_mean_141-150.mat']);
load([dpth,'output_h_mean_141-150.mat']);
load([dpth,'output_u_mean_141-150.mat']);
load([dpth,'output_v_mean_141-150.mat']);

h(1,:,:) = h(1,:,:) + eta;

% three layer ================================================
subplot(1,6,s)
hold on

for i = 2
    tmp_u = squeeze(u(i,2:end,:))*100;
    
    % reveal striations
    cutoff_length = 400;
    tmp_us = striations_reveal(tmp_u,re,cutoff_length);
    
    % zonal mean
    tmp_u_zm = nanmean(tmp_u(1:nx,1:ny/2),1);
    tmp_us_zm = nanmean(tmp_us(1:nx,1:ny/2),1);
    
    y1 = yy(1:ny/2);
    
    x1 = zeros(size(y1));
    l1 = plot(x1,y1,'-k','linewidth',2);   

    l2 = plot(tmp_us_zm,y1,'-k','linewidth',2);  
    
    
    x2 =  tmp_us_zm;x2(tmp_us_zm<0)=0;
    patch([x1,fliplr(x2)],[y1,fliplr(y1)],'r');
    x2 =  tmp_us_zm;x2(tmp_us_zm>0)=0;
    patch([x1,fliplr(x2)],[y1,fliplr(y1)],'b');

    le = legend([l2],{'us'});
    set(le,'location','best','edgecolor','none','fontsize',15); 
   
    title(['S',num2str(s),': u_{str} on layer ',num2str(i)],'fontweight','bold');
    xlabel('Zonal velocity (cm/s)');ylabel('y (km)');
    set(gca,'xtick',-.3:.1:.3,'xticklabel',{'-3','-2','-1','0','1','2','3'}...
        ,'ytick',1.e6:1.e6:5.e6,'yticklabel',{'1000','2000','3000','4000','5000'});
    set(gca,'fontsize',12,'box','off','tickdir','out','linewidth',1);
    box on; grid on;
    
end   
end


