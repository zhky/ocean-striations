%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare open ocean and marginal seas
% -------------------------------------------------------------------------
% - take Pacific and SCS as example
% - ECCO2, 2001-2010, 0.25*0.25
% -------------------------------------------------------------------------
% by zhangyu 20200107
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath(genpath('/home/dell/Striations/'))

%% initial parameters

% water depth level
lev = 1;

% resolution 
re = 0.25;

nc=netcdf('UVEL.1440x720x50.200101.nc');
lat = nc{'LATITUDE_T'}(:);
lon = nc{'LONGITUDE_T'}(:);
h = nc{'DEPTH_T'}(:);

%% process
% load data 
load SS_ECCO2_global_u_20012010mean.mat
u0 = utm*100; clear utm

% highpass filtering
us1 = zeros(size(u0));
for i = 1:length(h)
    i
    cutoff_length = 400;
    us1(i,:,:) = striations_reveal(squeeze(u0(i,:,:)),re,cutoff_length);
end

% highpass filtering
us2 = zeros(size(u0));
for i = 1:length(h)
    i
    cutoff_length = 200;
    us2(i,:,:) = striations_reveal(squeeze(u0(i,:,:)),re,cutoff_length);
end

%% striations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whole region
lon_west = 100;
lon_east = 205;  
lat_south = 5;
lat_north = 30;

% lines
lon_w1 = [105,135]+0.125;  
lon_e1 = [122,200]+0.125;  
lat_s1 = [  5, 10]+0.125;
lat_n1 = [ 24, 30]+0.125;  

figure('color','w','position',[100 100 1400 800])

for j = 1:length(lon_w1)
        
    A = find(lon_w1(j)<=lon & lon<=lon_e1(j));lon_box = lon(A);
    B = find(lat_s1(j)<=lat & lat<=lat_n1(j));lat_box = lat(B);
        
    % ------------------------------------------------------------------- %
    ax = subplot(3,2,1);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    vel_contours = [-40:10:-10,-9:-2,-1.5:0.5:1.5,2:9,10:10:40];

    m_contourf(lon_box,lat_box,squeeze(u0(1,B,A)),vel_contours,'linestyle','none'); 
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');
    
    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none','tickdir','out');  
    title('(a) North Pacific: u_{0}','fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

    % ------------------------------------------------------------------- %
    ax = subplot(3,2,3);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    % vel_contours = [-30:5:-10,-9:-2,-1.5:0.5:1.5,2:9,10:5:30];

    m_contourf(lon_box,lat_box,squeeze(us1(1,B,A)),vel_contours,'linestyle','none');
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');
%     co = colorbar('fontsize',15,'position',[0.88 0.16 0.02 0.7]); 
%     title(co,'cm s^-^1');   
%     set(co, 'tickdir', 'out');  

    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none','tickdir','out');  
    title(['(c) North Pacific: u_{s} on  {\lambda} = 400 km'],'fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

    % ------------------------------------------------------------------- %
    ax = subplot(3,2,5);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    % vel_contours = [-30:5:-10,-9:-2,-1.5:0.5:1.5,2:9,10:5:30];

    m_contourf(lon_box,lat_box,squeeze(us2(1,B,A)),vel_contours,'linestyle','none');
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');
    
    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none','tickdir','out');  
    title(['(e) North Pacific: u_{s} on  {\lambda} = 200 km'],'fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

end


% whole region
lon_west = 255;
lon_east = 345;  
lat_south = 13;
lat_north = 35;

% lines
lon_w1 = [260,295]+0.125;  
lon_e1 = [280,340]+0.125;  
lat_s1 = [ 17, 17]+0.125;
lat_n1 = [ 30, 35]+0.125;  

for j = 1:length(lon_w1)
        
    A = find(lon_w1(j)<=lon & lon<=lon_e1(j));lon_box = lon(A);
    B = find(lat_s1(j)<=lat & lat<=lat_n1(j));lat_box = lat(B);
        
    % ----------------------------------------------------------------- %
    ax = subplot(322);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    vel_contours = [-40:10:-10,-9:-2,-1.5:0.5:1.5,2:9,10:10:40];

    m_contourf(lon_box,lat_box,squeeze(u0(1,B,A)),vel_contours,'linestyle','none'); 
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');

    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none','tickdir','out');  
    title('(b) North Atlantic: u_{0}','fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

    % ----------------------------------------------------------------- %
    ax = subplot(324);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    % vel_contours = [-30:5:-10,-9:-2,-1.5:0.5:1.5,2:9,10:5:30];

    m_contourf(lon_box,lat_box,squeeze(us1(1,B,A)),vel_contours,'linestyle','none');
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');
    co = colorbar('fontsize',15,'position',[0.94 0.16 0.015 0.7]); 
    title(co,'cm s^-^1');   
    set(co, 'tickdir', 'out');   

    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none','tickdir','out');  
    title(['(d) North Atlantic: u_{str} on  {\lambda} = 400 km'],'fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

    % ----------------------------------------------------------------- %
    ax = subplot(326);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    % vel_contours = [-30:5:-10,-9:-2,-1.5:0.5:1.5,2:9,10:5:30];

    m_contourf(lon_box,lat_box,squeeze(us2(1,B,A)),vel_contours,'linestyle','none');
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');

    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none','tickdir','out');  
    title(['(f) North Atlantic: u_{str} on  {\lambda} = 200 km'],'fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

end

%% PSD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% whole region
lon_west = 90;
lon_east = 260;  
lat_south = 0;
lat_north = 30;

% lines
lon_w1 = [102,135]+0.125;  
lon_e1 = [120,235]+0.125;  
lat_s1 = [  0, 6]+0.125;
lat_n1 = [ 24, 30]+0.125;  

figure('color','w','position',[100 100 1600 800])

for j = 1:length(lon_w1)
        
    A = find(lon_w1(j)<=lon & lon<=lon_e1(j));lon_box = lon(A);
    B = find(lat_s1(j)<=lat & lat<=lat_n1(j));lat_box = lat(B);
        
    % ----------------------------------------------------------------- %
    img_1 = squeeze(u0(1,B,A));
    % remove mean and trend
    img_1 = detrend(img_1 -  nanmean(nanmean(img_1)));

    % psd
    s = re*110; 
    [PSD_1, k_1, l_1] = psd_period_2d(img_1,'hann', [size(img_1,1) size(img_1,2)], [1/s 1/s]);

    % one side PSD
    % transform two-side to one-side PSD
    k_11 = k_1(length(k_1)/2+1:end)*2;
    l_11 = l_1(length(l_1)/2+1:end)*2;
    PSD_11 = PSD_1(length(k_1)/2+1:end,length(l_1)/2+1:end)*2;

    subplot(3,4,j)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
    caxis([0 2000])
    cmocean('matter',20);colorbar;
    if j == 1
        title({'SCS', '(a1) PSD of u_0'},'fontsize',10,'fontweight','bold');
    elseif j == 2
        title({' SNP','(b1) PSD of u_0'},'fontsize',10,'fontweight','bold');
    end
    if j==1;    ylabel('{\lambda}_y (km)','fontweight','bold');end
    set(gca,'tickdir','out','fontsize',12);

    % ----------------------------------------------------------------- %
    img_1 = squeeze(us1(1,B,A));
    % remove mean and trend
    img_1 = detrend(img_1 -  nanmean(nanmean(img_1)));

    % psd
    s = re*110; 
    [PSD_1, k_1, l_1] = psd_period_2d(img_1,'hann', [size(img_1,1) size(img_1,2)], [1/s 1/s]);
    
    % one side PSD
    % transform two-side to one-side PSD
    k_11 = k_1(length(k_1)/2+1:end)*2;
    l_11 = l_1(length(l_1)/2+1:end)*2;
    PSD_11 = PSD_1(length(k_1)/2+1:end,length(l_1)/2+1:end)*2;

    subplot(3,4,j+4)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
    caxis([0 500])
    cmocean('matter',20);colorbar;
    if j == 1
        title('(a2) PSD of u_{s} ({\lambda}=400 km)','fontsize',10,'fontweight','bold');
    elseif j == 2
        title('(b2) PSD of u_{s} ({\lambda}=400 km)','fontsize',10,'fontweight','bold');
    end
    if j==1;    ylabel('{\lambda}_y (km)','fontweight','bold');end
    set(gca,'tickdir','out','fontsize',12);

    % ----------------------------------------------------------------- %
    img_1 = squeeze(us2(1,B,A));
    % remove mean and trend
    img_1 = detrend(img_1 -  nanmean(nanmean(img_1)));

    % psd
    s = re*110; 
    [PSD_1, k_1, l_1] = psd_period_2d(img_1,'hann', [size(img_1,1) size(img_1,2)], [1/s 1/s]);
    
    % one side PSD
    % transform two-side to one-side PSD
    k_11 = k_1(length(k_1)/2+1:end)*2;
    l_11 = l_1(length(l_1)/2+1:end)*2;
    PSD_11 = PSD_1(length(k_1)/2+1:end,length(l_1)/2+1:end)*2;

    subplot(3,4,j+8)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
    caxis([0 20])
    cmocean('matter',20);colorbar;
    if j == 1
        title('(a3) PSD of u_{s} ({\lambda}=200 km)','fontsize',10,'fontweight','bold');
    elseif j == 2
        title('(b3) PSD of u_{s} ({\lambda}=200 km)','fontsize',10,'fontweight','bold');
    end
    xlabel('{\lambda}_x (km)','fontweight','bold');
    if j==1;    ylabel('{\lambda}_y (km)','fontweight','bold');end
    set(gca,'tickdir','out','fontsize',12);
 
end

% whole region
lon_west = 250;
lon_east = 360;  
lat_south = 15;
lat_north = 35;

% lines
lon_w1 = [260,295]+0.125;  
lon_e1 = [280,340]+0.125;  
lat_s1 = [ 16, 15]+0.125;
lat_n1 = [ 30, 35]+0.125;  

for j = 1:length(lon_w1)
        
    A = find(lon_w1(j)<=lon & lon<=lon_e1(j));lon_box = lon(A);
    B = find(lat_s1(j)<=lat & lat<=lat_n1(j));lat_box = lat(B);
        
    % ----------------------------------------------------------------- %
    img_1 = squeeze(u0(1,B,A));
    % remove mean and trend
    img_1 = detrend(img_1 -  nanmean(nanmean(img_1)));

    % psd
    s = re*110; 
    [PSD_1, k_1, l_1] = psd_period_2d(img_1,'hann', [size(img_1,1) size(img_1,2)], [1/s 1/s]);

    % one side PSD
    % transform two-side to one-side PSD
    k_11 = k_1(length(k_1)/2+1:end)*2;
    l_11 = l_1(length(l_1)/2+1:end)*2;
    PSD_11 = PSD_1(length(k_1)/2+1:end,length(l_1)/2+1:end)*2;

    subplot(3,4,j+2)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
%     caxis([0 1000])
    cmocean('matter',20);colorbar;
    if j == 1
        title({'GM','(c1) PSD of u_0'},'fontsize',10,'fontweight','bold');
    elseif j == 2
        title({'NA','(d1) PSD of u_0'},'fontsize',10,'fontweight','bold');
    end
    set(gca,'tickdir','out','fontsize',12);

    % ----------------------------------------------------------------- %
    img_1 = squeeze(us1(1,B,A));
    % remove mean and trend
    img_1 = detrend(img_1 -  nanmean(nanmean(img_1)));

    % psd
    s = re*110; 
    [PSD_1, k_1, l_1] = psd_period_2d(img_1,'hann', [size(img_1,1) size(img_1,2)], [1/s 1/s]);
    
    % one side PSD
    % transform two-side to one-side PSD
    k_11 = k_1(length(k_1)/2+1:end)*2;
    l_11 = l_1(length(l_1)/2+1:end)*2;
    PSD_11 = PSD_1(length(k_1)/2+1:end,length(l_1)/2+1:end)*2;

    subplot(3,4,j+6)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
%     caxis([0 200])
    cmocean('matter',20);colorbar;
    if j == 1
        title('(c2) PSD of u_{str} ({\lambda}=400 km)','fontsize',10,'fontweight','bold');
    elseif j == 2
        title('(d2) PSD of u_{str} ({\lambda}=400 km)','fontsize',10,'fontweight','bold');
    end
    set(gca,'tickdir','out','fontsize',12);

    % ----------------------------------------------------------------- %
    img_1 = squeeze(us2(1,B,A));
    % remove mean and trend
    img_1 = detrend(img_1 -  nanmean(nanmean(img_1)));

    % psd
    s = re*110; 
    [PSD_1, k_1, l_1] = psd_period_2d(img_1,'hann', [size(img_1,1) size(img_1,2)], [1/s 1/s]);
    
    % one side PSD
    % transform two-side to one-side PSD
    k_11 = k_1(length(k_1)/2+1:end)*2;
    l_11 = l_1(length(l_1)/2+1:end)*2;
    PSD_11 = PSD_1(length(k_1)/2+1:end,length(l_1)/2+1:end)*2;

    subplot(3,4,j+10)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
%     caxis([0 10])
    cmocean('matter',20);colorbar;
    if j == 1
        title('(c3) PSD of u_{str} ({\lambda}=200 km)','fontsize',10,'fontweight','bold');
    elseif j == 2
        title('(d3) PSD of u_{str} ({\lambda}=200 km)','fontsize',10,'fontweight','bold');
    end
    xlabel('{\lambda}_x (km)','fontweight','bold');
    set(gca,'tickdir','out','fontsize',12);
 
end
