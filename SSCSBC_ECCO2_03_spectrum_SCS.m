%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare open ocean and marginal seas
% -------------------------------------------------------------------------
% - take Pacific and SCS as example
% - ECCO2, 2001-2010, 0.25*0.25
% -------------------------------------------------------------------------
% by zhangyu 20200107
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initial parameters
clear
addpath(genpath('/home/dell/01_PROJECTS/'));

% water depth level
lev = 1;

% resolution 
re = 0.25;

% data path
str = '/home/dell/data/02_MODEL/ECCO2/UVEL';

%number of files 
file_stru=dir(['',str,'']);
num=size(file_stru,1);
filename=strcat(['',str,'','/',file_stru(3).name]);
nc=netcdf(filename,'nowrite');
lat=nc{'LATITUDE_T'}(:);
lon=nc{'LONGITUDE_T'}(:);       
h=nc{'DEPTH_T'}(:);

close(nc)

%% process
% load data 
load SSCS_ECCO2_global_u_20012010mean.mat
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

figure('color','w','position',[50 50 1200 900]);

for j = 1:length(lon_w1)
        
    A = find(lon_w1(j)<=lon & lon<=lon_e1(j));lon_box = lon(A);
    B = find(lat_s1(j)<=lat & lat<=lat_n1(j));lat_box = lat(B);
        
    % ----------------------------------------------------------------- %
    ax = subplot(311);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    vel_contours = [-50,-30:5:-10,-9:-2,-1.5:0.5:1.5,2:9,10:5:30,50];

    m_contourf(lon_box,lat_box,squeeze(u0(1,B,A)),vel_contours,'linestyle','none'); 
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');
    co = colorbar('fontsize',15); title(co,'cm s^-^1');   

    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none');  
    title('(a)North Pacific: u_{0}','fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

    % ----------------------------------------------------------------- %
    ax = subplot(312);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    % vel_contours = [-30:5:-10,-9:-2,-1.5:0.5:1.5,2:9,10:5:30];

    m_contourf(lon_box,lat_box,squeeze(us1(1,B,A)),vel_contours,'linestyle','none');
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');
    co = colorbar('fontsize',15); title(co,'cm s^-^1');   

    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none');  
    title(['(b)North Pacific: u_{str} on  {\lambda} = 400 km'],'fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

    % ----------------------------------------------------------------- %
    ax = subplot(313);
    hold on
    m_proj('robinson','lon',[lon_west lon_east],'lat',[lat_south lat_north]);

    % vel_contours = [-30:5:-10,-9:-2,-1.5:0.5:1.5,2:9,10:5:30];

    m_contourf(lon_box,lat_box,squeeze(us2(1,B,A)),vel_contours,'linestyle','none');
    ax = customcmap_axs(ax,vel_contours,create_colors(length(vel_contours)-1,[0 0 0.2; 0 0 1; 1 1 1; 1 0 0; 0.2 0 0]),...
        [min(vel_contours) max(vel_contours)]);
    colorbar('off');
    co = colorbar('fontsize',15); title(co,'cm s^-^1');   

    m_coast('patch',[.7 .7 .7],'edgecolor','none','linewidth',1);
    m_grid('fontsize',12,'color','k','linestyle','none');  
    title(['(c)North Pacific: u_{str} on  {\lambda} = 200 km'],'fontsize',15,'fontweight','bold');

    % add boxes
    for i = 1:length(lon_w1)
        m_line([lon_w1(i),lon_e1(i),lon_e1(i),lon_w1(i),lon_w1(i)],...
            [lat_s1(i),lat_s1(i),lat_n1(i),lat_n1(i),lat_s1(i)],'color','k','linewidth',2);    
    end

end

%% PSD
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

figure('color','w','position',[50 50 1200 900]);

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

    subplot(3,2,j)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
%     caxis([0 1000])
    cmocean('matter',20);colorbar;
    if j == 1
        title('2D psd of u_0 in SCS','fontsize',12,'fontweight','bold');
    elseif j == 2
        title('2D psd of u_0 in SNP','fontsize',12,'fontweight','bold');
    end
    ylabel('L in meridional (km)','fontweight','bold');
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

    subplot(3,2,j+2)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
%     caxis([0 200])
    cmocean('matter',20);colorbar;
    if j == 1
        title('2D psd of u_{str} on {\lambda} = 400 km in SCS','fontsize',12,'fontweight','bold');
    elseif j == 2
        title('2D psd of u_{str} on {\lambda} = 400 km in SNP','fontsize',12,'fontweight','bold');
    end
    ylabel('L in meridional (km)','fontweight','bold');
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

    subplot(3,2,j+4)
    contourf(1./l_11,1./k_11,PSD_11,100,'linestyle','none');
%     caxis([0 10])
    cmocean('matter',20);colorbar;
    if j == 1
        title('2D psd of u_{str} on {\lambda} = 200 km in SCS','fontsize',12,'fontweight','bold');
    elseif j == 2
        title('2D psd of u_{str} on {\lambda} = 200 km in SNP','fontsize',12,'fontweight','bold');
    end
    xlabel('L in zonal (km)','fontweight','bold');
    ylabel('L in meridional (km)','fontweight','bold');
    set(gca,'tickdir','out','fontsize',12);
 
end
