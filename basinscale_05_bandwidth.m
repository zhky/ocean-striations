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

pth0 = ['/home/dell/Striations/striations_3layers_basinscale/'];

R = 6371*1.e3;      % radius of the earth, unit: m

% basinscale_3layer_s0
xlen0 = [2880, 3840, 5760, 7680];
ylen0 = [3840, 3840, 3840, 3840];
nx0 =   [ 97, 129, 193, 257];
ny0 =   [129, 129, 129, 129];

% basinscale_3layer_re15km
xlen1 = [ 960, 1920, 2880, 3840, 4800, 5760];
ylen1 = [3840, 3840, 3840, 3840, 3840, 5760];
nx1 = [ 64, 128, 192, 256, 320, 384];
ny1 = [256, 256, 256, 256, 256, 256];


figure('color','w','position',[100 100 1200 1000])
%% four scenaros 
for s = 1:6
xlen = xlen1(s)*1.e3;
ylen = ylen1(s)*1.e3;
nx = nx1(s);
ny = ny1(s);
dx = xlen/nx;
dy = ylen/ny;
re = dx/1.e3/100;

dpth = [pth0,'s1', num2str(s),'/'];

load([dpth,'output_eta_mean_141-150.mat']);
load([dpth,'output_h_mean_141-150.mat']);
load([dpth,'output_u_mean_141-150.mat']);
load([dpth,'output_v_mean_141-150.mat']);

h(1,:,:) = h(1,:,:) + eta;

% u_0 and u_str ================================================
%figure('color','w','position',[100 100 1500 800])

for i = 1:3
    % u_0 -----------------------------------------
    tmp_u = squeeze(u(i,2:end,:));%size(tmp_u)
    tmp_v = squeeze(v(i,:,2:end));%size(tmp_v)
    [X,Y] = meshgrid(dx:dx:xlen,dy:dy:ylen);
    
    % u_str -------------------------------------
    cutoff_length = 400;
    tmp_us = striations_reveal(tmp_u,re,cutoff_length);

    % eastward jet detecting
    % jet parameters
    re_x = dx;
    re_y = dy;
    jet_distinguish = 5*re_y;  % separate two jets whose distance is large than jet_distinguish
    jet_length = 1.5/re_x;       % save the jets whose length is longer than jet_length

    x1 = dx:dx:xlen;
    y1 = dy:dy:ylen;
    data1 = tmp_us*10;data1(data1<0.005)=0;
    [wholejet_value1,wholejet_x1,wholejet_y1,jetaxis_value1,jetaxis_x1,jetaxis_y1] = ...
        find_wholejet_02(data1',x1,y1,jet_distinguish,jet_length);

    % westward jet detecting
    data2 = tmp_us*10;data2(data2>-0.005)=0;
    [wholejet_value2,wholejet_x2,wholejet_y2,jetaxis_value2,jetaxis_x2,jetaxis_y2] = ...
        find_wholejet_02(data2',x1,y1,jet_distinguish,jet_length);

    % all
    jetaxis_x=[jetaxis_x1;jetaxis_x2];
    jetaxis_y=[jetaxis_y1;jetaxis_y2];
    wholejet_value=[wholejet_value1;wholejet_value2];
    wholejet_x=[wholejet_x1;wholejet_x2];
    wholejet_y=[wholejet_y1;wholejet_y2];

    % bandwidth
    bandwidth = zeros(size(wholejet_value,1),size(wholejet_value,2))*nan;
    for bi = 1:size(wholejet_value,1)
        for bj = 1:size(wholejet_value,2)
            tmp_y = squeeze(wholejet_y(bi,bj,:));
            tmp_v = squeeze(wholejet_value(bi,bj,:));
            tmp_y(tmp_v == 0) = [];
            if length(tmp_y) >=2
                 bandwidth(bi,bj) = tmp_y(end) - tmp_y(1);
            end
        end
    end

    % plot bandwidth statistic ----------------------------------------
    %subplot(3,4,(s-3)+(i-1)*4)
    subplot(3,2,i*2-1)
    hold on
    
    bandwidth(bandwidth<=0) = nan;
    data = bandwidth(:)/1000;
    data(isnan(data)) = [];

    % 检验是否符合正态分布
    % [h, p] = kstest(data);
    
    % 计算平均值和标准差
    mu = nanmean(data);
    sigma = nanstd(data);
    data(data > mu+3*sigma) = [];
    eval(['mu_s_l',num2str(i),'(s) = mu;']);
    
    % 绘频率分布曲线
    %h1 = histogram(data,'Normalization','pdf');
    
    counts = histcounts(data, 10); % 默认10个bins
    BinEdges = linspace(min(data), max(data), 11); % 获取bin的边界
    L = length(BinEdges);
    hx1 = zeros(1,L-1);
    hy1 = zeros(1,L-1);
    for j = 1:L-1
	hx1(j) = (BinEdges(j) + BinEdges(j+1))/2;
    end
    hy1 = smooth(counts,5);

    eval(['l',num2str(s),'= plot(hx1,hy1,''linewidth'',2);']);
    %plot(mu,hy1(hx1==mu),'.k','markersize',30);
    %plot([mu mu],[0 nanmax(hy1)],'--k');
    %axis([0 300 0 5.e-3]);
    hold on; 

    ylabel('Counts','Fontsize',15);
    if i == 1
    	title(['(a) layer ',num2str(i)],'fontweight','bold','fontsize',15);
    elseif i == 2
	title(['(c) layer ',num2str(i)],'fontweight','bold','fontsize',15);
    else
    	title(['(e) layer ',num2str(i)],'fontweight','bold','fontsize',15);
    end
    box on; grid on;axis tight;
    

end
end

xlabel('Statistical striations width (km)','Fontsize',15);
le = legend([l1,l2,l3,l4,l5,l6],'s1','s2','s3','s4','s5','s6');
set(le,'location','best','fontsize',15,'edgecolor','none');

% plot bandwidth rise of basinwidth
for i = 1:3
subplot(3,2,i*2)
eval(['plot(mu_s_l',num2str(i),',''-ok'',''linewidth'',2);']);
if i == 1
    	title(['(b) layer ',num2str(i)],'fontweight','bold','fontsize',15);
    elseif i == 2
	title(['(d) layer ',num2str(i)],'fontweight','bold','fontsize',15);
    else
    	title(['(f) layer ',num2str(i)],'fontweight','bold','fontsize',15);
end
ylabel('Mean width (km)','fontsize',15);
box on; grid on;axis tight;

set(gca,'xtick',1:6,'xticklabel',{'s1','s2','s3','s4','s5','s6'});
end
xlabel('Scenarios','Fontsize',15);

mu_s_l1
mu_s_l2
mu_s_l3

