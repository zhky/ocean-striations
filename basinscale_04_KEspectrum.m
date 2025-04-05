%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to verify the optimal layer view of striations by using idea isopycnal
% model "Aronnax".
% check the flow field after highpass filtering
%  ------------------------------------------------------------------------
% by zhangyu 2020526
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
xlen1 = [ 960, 1920, 2880, 3840, 4800, 5760];
ylen1 = [3840, 3840, 3840, 3840, 3840, 3840];
nx1 = [ 64, 128, 192, 256, 320, 384];
ny1 = [256, 256, 256, 256, 256, 256];

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
figure('color','w','position',[100 100 1500 800])

for s = 1:6  % loop of scenaros
xlen = xlen1(s)*1.e3;
ylen = ylen1(s)*1.e3;
nx = nx1(s);
ny = ny1(s);
dx = xlen/nx;
dy = ylen/ny;
re = dx/1.e3/100;

dpth = [pth0,'s1',num2str(s),'/'];

load([dpth,'output_eta_mean.mat']);
load([dpth,'output_h_mean.mat']);
load([dpth,'output_u_mean.mat']);
load([dpth,'output_v_mean.mat']);

h(1,:,:) = h(1,:,:) + eta;

% u_0 and u_str ================================================
%figure('color','w','position',[100 100 1500 800])

for i = 1:3  % loop of layer
    % u_0
    tmp_u0 = squeeze(u(i,2:end,:));
    tmp_v0 = squeeze(v(i,:,2:end));
    [X,Y] = meshgrid(dx:dx:xlen,dy:dy:ylen);
    
    % u_str -------------------------------------
    cutoff_length = 400;
    tmp_us = striations_reveal(tmp_u0,re,cutoff_length);
      
    % remove mean and trend
    usa = detrend(tmp_us -  nanmean(nanmean(tmp_us)));
    u0a = detrend(tmp_u0 -  nanmean(nanmean(tmp_u0)));

    % psd --------------------------------
    sp = re*110;
    [PSD_1, k_1, l_1] = psd_period_2d...
        (usa','hann', [size(usa,1) size(usa,2)], [1/sp 1/sp]);
    [PSD_0, k_0, l_0] = psd_period_2d...
        (u0a','hann', [size(u0a,1) size(u0a,2)], [1/sp 1/sp]);

    % one side PSD
    % transform two-side to one-side PSD
    k_11 = k_1(length(k_1)/2+1:end)*2;
    l_11 = l_1(length(l_1)/2+1:end)*2;
    PSD_11 = PSD_1(length(k_1)/2+1:end,length(l_1)/2+1:end)*2;

    k_00 = k_0(length(k_0)/2+1:end)*2;
    l_00 = l_0(length(l_0)/2+1:end)*2;
    PSD_00 = PSD_0(length(k_0)/2+1:end,length(l_0)/2+1:end)*2;

    %{
    subplot(length(targ_h),3,1+3*(h-1))
    contourf(1./l_00,1./k_00,PSD_00,100,'linestyle','none');
    %caxis([0 200])
    cmocean('matter',20);colorbar;
    title(['ECCO2_-U0: ',num2str(targ_h(h)),'m'],'fontsize',15);
    ylabel('L in meridional (km)','fontweight','bold');
    set(gca,'tickdir','out','fontsize',12);
    %}

    subplot(3,6,(s)+(i-1)*6)
    cf = [0:1.e-7:1.e-6,2e-6];
    contourf(1./l_11,1./k_11,PSD_11,cf,'linestyle','none');
    caxis([0 1.e-6])
    cmocean('matter',20);

    if s==6 & i==2
    cb = colorbar('position',[0.92 0.15 0.02 0.7]);
    set(get(cb,'title'),'string',' ');
    set(cb,'tickdir','out')
    end

    title(['S',num2str(s),': layer',num2str(i)],'fontsize',15);
    if i==3;xlabel('L_x (km)','fontweight','bold');end
    if s==1;ylabel('L_y (km)','fontweight','bold');end
    set(gca,'tickdir','out','fontsize',12);

    end % for i
end % for s


