%% This code will plot the variance field of temperature and salinity from
% the runs. 
% NOTE OF CAUTION: If there is not enough data, be careful in interpreting
% these results
% COMMENT THESE LINES IF YOUR USING THIS FILE WITH THE SHELL SCRIPT
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%% File directory name (Fill the location and name of the file seperately)
file_dir = pwd;
list_temp = dir(fullfile(file_dir,'T.0*.data'));
list_salt = dir(fullfile(file_dir,'S.0*.data'));
list_U = dir(fullfile(file_dir,'U*.data'));
list_V = dir(fullfile(file_dir,'V*.data'));
list_W = dir(fullfile(file_dir,'W*.data'));
list_surf = dir(fullfile(file_dir,'surf*.data'));
topofile = [file_dir filesep 'shelficeTopo.Lin.bin'];
delta_t = 4;
day2s = 24 * 3600;
day2hr = 24;
time_start = 22; % In days
time_end = 23;% In days
%% Co-ordinate files
x_max = 4000;
Y = cumsum(squeeze(rdmds('DYC')));
Z = -cumsum(squeeze(rdmds('DRF')));
Z_start = -100; % The minimum depth of ice-shelf (approx.)
Z_index = find(Z>Z_start,1,'last');
Z(1:Z_index) = [];
[Ymesh,Zmesh] = meshgrid(Y,Z);

%% Skipping grid points in both direction
hor_skip = 4;
ver_skip = 1;

%% Timestamps of 2D fields
for i=1:length(list_salt)
    time_2D(i) = str2num(list_salt(i).name(end-10:end-5)) * delta_t;
end

%% Timestamps of the interface properties
for i=1:length(list_surf)
    time_interface(i) = str2num(list_surf(i).name(end-10:end-5)) * delta_t;
end

%% Extracting the interface properties
for i=1:length(list_surf)
    Data_2D = squeeze(rdmds(list_surf(i).name(1:end-5)));
    Um(i,:) = Data_2D(:,end);
end

%% Interpolating in time
for i=1:size(Um,2)
    Um_i(:,i) = interp1(time_interface,Um(:,i),time_2D);
end

Ymesh = Ymesh(1:hor_skip:end,1:ver_skip:end);
Zmesh = Zmesh(1:hor_skip:end,1:ver_skip:end);

%% Extracting the time stamps which are requested
time_index = find(time_2D>=time_start*day2s&time_2D<=time_end*day2s);

%% Finding the mean of the temperature and salinity field
mean.T = zeros(size(Ymesh'));
mean.S = zeros(size(Zmesh'));
mean.melt_surf = zeros(size(Ymesh(1,:)))';
% mean.V = zeros(size(Ymesh'));
% mean.W = zeros(size(Zmesh'));
count = 1;
for i=time_index
    
    T = squeeze(rdmds(list_temp(i).name(1:end-5)));
    S = squeeze(rdmds(list_salt(i).name(1:end-5)));
    S(S==0) = NaN;
    T(isnan(S)) = NaN;
    
    T(:,1:Z_index) = [];
    S(:,1:Z_index) = [];
    
    T = T(1:ver_skip:end,1:hor_skip:end);
    S = S(1:ver_skip:end,1:hor_skip:end);
    
    for j=1:size(S,1)
        try
            Smm(j) = S(j,find(~isnan(S(j,:)),1,'first'));
            Tmm(j) = T(j,find(~isnan(T(j,:)),1,'first'));
        catch
            Smm(j) = NaN;
            Tmm(j) = NaN;
        end
    end
    
    mean.T = mean.T + T;
    mean.S = mean.S + S;
%     mean.Vel = mean.Vel + Vel;
%     [Tb,Sb,melt(count,:)] = meltrate_gamm_dep_new(Tmm',Smm',Umm',topofile,0);
    [Tbb,Sbb,melt1(count,:)] = meltrate_gamm_dep_new(Tmm',Smm',Um_i(i,:)',topofile,0);
    %mean.melt = mean.melt + melt;
    %mean.melt_surf = mean.melt_surf + melt1;
    count = count + 1;
end

mean.T = mean.T/(length(time_index));
mean.S = mean.S/(length(time_index));
%mean.Vel = mean.Vel/(length(list_temp) - 1);
% mean.melt = mean.melt/(length(list_temp) - 1);
% mean.melt_surf = mean.melt_surf/(length(list_temp) - 1);
%mean.melt = nanmean(melt);
mean.melt_surf = nanmean(melt1);

% clear Um;
% for j=1:size(mean.S,1)
%     try
%         Sm(j) = mean.S(j,find(~isnan(mean.S(j,:)),1,'first'));
%         Tm(j) = mean.T(j,find(~isnan(mean.T(j,:)),1,'first'));
%         Um(j) = mean.Vel(j,find(~isnan(mean.Vel(j,:)),1,'first'));
%     catch
%         Sm(j) = NaN;
%         Tm(j) = NaN;
%         Um(j) = NaN;
%     end
% end
% [Tb,Sb,wb] = meltrate_gamm_dep_new(Tm',Sm',Um',topofile,0);

%% Plotting the mean fields as a figure
figure('units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(3,1,1);
contourf(ax1,Ymesh,Zmesh,mean.T');
colorbar;
colormap(ax1,'winter');
xlabel(ax1,'Distance from GL (m)');
zlabel(ax1,'Depth (m)');
title(ax1,'Temperature field');

ax2 = subplot(3,1,2);
contourf(ax2,Ymesh,Zmesh,mean.S');
colorbar;
colormap(ax2,'jet');
xlabel(ax2,'Distance from GL (m)');
zlabel(ax2,'Depth (m)');
title(ax2,'Salinity field');

% ax3 = subplot(5,1,3);
% ax3_p = ax3.InnerPosition;
% ax2_p = ax2.InnerPosition;
% set(ax3,'InnerPosition',[ax3_p(1) ax3_p(2) ax2_p(3) ax3_p(4)]);
% %wb(Ymesh(1,:)>=x_max) = 0;
% plot(ax3,Ymesh,wb);
% %set(ax3,'XTick',ax2,XTick);
% xlabel(ax3,'Distance from GL (m)');
% ylabel(ax3,'$\dot{m} (m/yr)$','interpreter','latex');
% %xlim(ax3,[0 x_max]);
% title(ax3,'Melt rate');
% 
% ax4 = subplot(5,1,4);
% ax4_p = ax4.InnerPosition;
% ax2_p = ax2.InnerPosition;
% set(ax4,'InnerPosition',[ax4_p(1) ax4_p(2) ax2_p(3) ax4_p(4)]);
% %wb(Ymesh(1,:)>=x_max) = 0;
% plot(ax4,Ymesh,mean.melt);
% %set(ax3,'XTick',ax2,XTick);
% xlabel(ax4,'Distance from GL (m)');
% ylabel(ax4,'$\dot{m} (m/yr)$','interpreter','latex');
% %xlim(ax3,[0 x_max]);
% title(ax4,'Melt rate (averaged)');

ax5 = subplot(3,1,3);
ax5_p = ax5.InnerPosition;
ax2_p = ax2.InnerPosition;
set(ax5,'InnerPosition',[ax5_p(1) ax5_p(2) ax2_p(3) ax5_p(4)]);
%wb(Ymesh(1,:)>=x_max) = 0;
plot(ax5,Ymesh,mean.melt_surf);
%set(ax3,'XTick',ax2,XTick);
xlabel(ax5,'Distance from GL (m)');
ylabel(ax5,'$\dot{m} (m/yr)$','interpreter','latex');
%xlim(ax3,[0 x_max]);
title(ax5,'Melt rate');

suptitle(sprintf('Mean fields (Day %d - Day %d)',time_start,time_end));
print(gcf,'-dpng','-r300',fullfile(pwd,sprintf('Mean_fields_%d_%d',time_start,time_end)));