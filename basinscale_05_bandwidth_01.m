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
xlen1 = [ 960, 1920, 2880, 3840, 4800];
ylen1 = [3840, 3840, 3840, 3840, 3840];
nx1 = [ 64, 128, 192, 256, 320];
ny1 = [256, 256, 256, 256, 256];


figure('color','w','position',[100 100 1200 900])
%% four scenaros 
for s = 1:5
xlen = xlen1(s)*1.e3;
ylen = ylen1(s)*1.e3;
nx = nx1(s);
ny = ny1(s);
dx = xlen/nx;
dy = ylen/ny;
re = dx/1.e3/100;

dpth = [pth0,'s1', num2str(s),'/'];

load([dpth,'output_eta_mean.mat']);
load([dpth,'output_h_mean.mat']);
load([dpth,'output_u_mean.mat']);
load([dpth,'output_v_mean.mat']);

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

    [h, p] = kstest(data);
    
    % 计算平均值和标准差
    mu = nanmean(data);
    sigma = nanstd(data);
    data(data > mu+3*sigma) = [];

    eval(['mu_s_l',num2str(i),'(s) = mu;']);
    
    % 绘制标准正态分布曲线
    x = linspace(0, 500, 100);
    y = normpdf(x, mu, sigma);
    %histogram(data,30,'Normalization','pdf');
    eval(['l',num2str(s),'= plot(x,y,''linewidth'',2);']);
    plot(mu,nanmax(y),'.k','markersize',30);
    plot([mu mu],[0 nanmax(y)],'--k');
    hold on; 

    %{
    % 在x轴上标示出距离平均值位置的1倍、2倍和3倍标准差值位置
    line([mu-sigma, mu-sigma], [0, normpdf(mu-sigma, mu, sigma)], ...
         'LineStyle', '--', 'LineWidth', 2);
    line([mu+sigma, mu+sigma], [0, normpdf(mu+sigma, mu, sigma)], ...
         'LineStyle', '--', 'LineWidth', 2); 
    legend([l1,l2],'标准正态分布', '1倍标准差');
    %}

    %xlabel('Normal Distribution of Striations Bandwidth (km)','Fontsize',15);
    ylabel('Probability Denstiy','Fontsize',15);
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

xlabel('Striations width (km)','Fontsize',15);
le = legend([l1,l2,l3,l4,l5],'s1','s2','s3','s4','s5');
set(le,'location','best','fontsize',15,'edgecolor','none');

% plot bandwidth rise of basinwidth
for i = 1:3
subplot(3,2,i*2)
eval(['plot(mu_s_l',num2str(i),',''-o'',''linewidth'',2);']);
if i == 1
    	title(['(b) layer ',num2str(i)],'fontweight','bold','fontsize',15);
    elseif i == 2
	title(['(d) layer ',num2str(i)],'fontweight','bold','fontsize',15);
    else
    	title(['(f) layer ',num2str(i)],'fontweight','bold','fontsize',15);
end
ylabel('Striations width (km)','fontsize',15);
box on; grid on;axis tight;

xlabel('Scenarios','Fontsize',15);
set(gca,'xtick',1:5,'xticklabel',{'s1','s2','s3','s4','s5'});
end
