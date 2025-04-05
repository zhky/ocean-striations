%%%%%%%%%%%%%%%%i%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% figure('color','w','position',[100 100 800 800])

%% four scenaros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s = 1
xlen = xlen1(s)*1.e3;
ylen = ylen1(s)*1.e3;
nx = nx1(s);
ny = ny1(s);
dx = xlen/nx;
dy = ylen/ny;
re = dx/1.e3/100;

dpth = [pth0,'/s1', num2str(s),'/'];

load([dpth,'output_eta_mean.mat']);
load([dpth,'output_h_mean.mat']);
load([dpth,'output_u_mean.mat']);
load([dpth,'output_v_mean.mat']);

h(1,:,:) = h(1,:,:) + eta;

% u_0 and u_str ================================================
%figure('color','w','position',[100 100 800 1000])

for k = 2
    % u_0 -----------------------------------------
    tmp_u = squeeze(u(k,2:end,:));%size(tmp_u)
    tmp_v = squeeze(v(k,:,2:end));%size(tmp_v)
    [X,Y] = meshgrid(dx:dx:xlen,dy:dy:ylen);
    
    % u_str -------------------------------------
    cutoff_length = 400;
    tmp_us = striations_reveal(tmp_u,re,cutoff_length);

    % eastward jet detecting
    % jet parameters
    re_x = dx;
    re_y = dy;
    jet_distinguish = 2*re_y;  % separate two jets whose distance is large than jet_distinguish
    jet_length = 1.5/re_x;       % save the jets whose length is longer than jet_length

    x1 = dx:dx:xlen;
    y1 = dy:dy:ylen;
    data1 = tmp_us*10;data1(data1<0.0005)=0;
    [wholejet_value1,wholejet_x1,wholejet_y1,jetaxis_value1,jetaxis_x1,jetaxis_y1] = ...
        find_wholejet_02(data1',x1,y1,jet_distinguish,jet_length);

    % westward jet detecting
    data2 = tmp_us*10;data2(data2>-0.0005)=0;
    [wholejet_value2,wholejet_x2,wholejet_y2,jetaxis_value2,jetaxis_x2,jetaxis_y2] = ...
        find_wholejet_02(data2',x1,y1,jet_distinguish,jet_length);

    % all
    jetaxis_x=[jetaxis_x1;jetaxis_x2];
    jetaxis_y=[jetaxis_y1;jetaxis_y2];
    wholejet_value=[wholejet_value1;wholejet_value2];
    wholejet_x=[wholejet_x1;wholejet_x2];
    wholejet_y=[wholejet_y1;wholejet_y2];

    % all
    jetaxis_x=[jetaxis_x1;jetaxis_x2];
    jetaxis_y=[jetaxis_y1;jetaxis_y2];
    wholejet_value=[wholejet_value1;wholejet_value2];
    wholejet_x=[wholejet_x1;wholejet_x2];
    wholejet_y=[wholejet_y1;wholejet_y2];
    
    % plot axes and bandwidth -------------------------------------------------------
    %subplot(1,3,s)
    figure('color','w','position',[100 100 800 600])
    hold on
    contourf(X,Y,tmp_us',100,'linestyle','none');
    cb = colorbar;
    set(get(cb,'title'),'string','m s^{-1}')

    colormap(create_colors(100-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]));
    caxis([-0.005 0.005]);

    %{
    title(['S',num2str(s),': u_{s} on layer ',num2str(i)],'fontweight','bold');
    xlabel('x (1e3 km)');ylabel('y (1e3 km)');
    set(gca,'xtick',1.e6:1.e6:8.e6,'xticklabel',{'1','2','3','4','5','6','7','8'}...
        ,'ytick',1.e6:1.e6:8.e6,'yticklabel',{'1','2','3','4','5','6','7','8'});
    set(gca,'fontsize',12,'box','off','tickdir','out','linewidth',1);
    %}

    disp('ploting jets band')
    k1 = 1;
    for i = 1:size(wholejet_value,1)
        for j = 1:k1:size(wholejet_value,2)
            tmp_x = squeeze(wholejet_x(i,j,:));
            tmp_y = squeeze(wholejet_y(i,j,:));
            tmp_v = squeeze(wholejet_value(i,j,:));
            tmp_x(tmp_v == 0) = [];
            tmp_y(tmp_v == 0) = [];
            line(tmp_x,tmp_y,'linestyle','-','color',[.5 .5 .5],'linewidth',1.5);
        end
    end

    disp('ploting jets axes')
    for i = 1:size(jetaxis_y,1)
        tmpx = jetaxis_x(i,:);tmpx(tmpx==0) = [];
        tmpy = jetaxis_y(i,:);tmpy(tmpy==0) = [];
        line(tmpx,tmpy,'linestyle','-','color','k','linewidth',1.5);
    end

    title(['S',num2str(s),': u_{s} on layer ',num2str(k)],'fontweight','bold');
    xlabel('x (km)');ylabel('y (km)');
    set(gca,'xtick',2.e5:2.e5:1.e6,'xticklabel',{'200','400','600','800','1000'}...
        ,'ytick',1.e6:1.e6:4.e6,'yticklabel',{'1000','2000','3000','4000'});
    set(gca,'fontsize',12,'box','on','tickdir','out','linewidth',2);


   end
end

