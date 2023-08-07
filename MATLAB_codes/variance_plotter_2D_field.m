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
topofile = [file_dir filesep 'shelficeTopo.Lin.bin'];

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

Ymesh = Ymesh(1:hor_skip:end,1:ver_skip:end);
Zmesh = Zmesh(1:hor_skip:end,1:ver_skip:end);

%% Finding the mean of the temperature and salinity field
mean.T = zeros(size(Ymesh'));
mean.S = zeros(size(Zmesh'));
mean.Vel = zeros(size(Ymesh'));
% mean.V = zeros(size(Ymesh'));
% mean.W = zeros(size(Zmesh'));

for i=1:length(list_temp)
    
    T = squeeze(rdmds(list_temp(i).name(1:end-5)));
    S = squeeze(rdmds(list_salt(i).name(1:end-5)));
    U = squeeze(rdmds(list_U(i).name(1:end-5)));
    V = squeeze(rdmds(list_V(i).name(1:end-5)));
    W = squeeze(rdmds(list_W(i).name(1:end-5)));
    S(S==0) = NaN;
    T(isnan(S)) = NaN;
    U(isnan(S)) = NaN;
    V(isnan(S)) = NaN;
    W(isnan(S)) = NaN;
    Vel = sqrt(U.^2 + V.^2 + W.^2);
    
    T(:,1:Z_index) = [];
    S(:,1:Z_index) = [];
    Vel(:,1:Z_index) = [];
    
    T = T(1:ver_skip:end,1:hor_skip:end);
    S = S(1:ver_skip:end,1:hor_skip:end);
    Vel = Vel(1:ver_skip:end,1:hor_skip:end);
    
    mean.T = mean.T + T;
    mean.S = mean.S + S;
    mean.Vel = mean.Vel + Vel;
end

mean.T = mean.T/(length(list_temp) - 1);
mean.S = mean.S/(length(list_salt) - 1);
mean.Vel = mean.Vel/(length(list_temp) - 1);

for j=1:size(mean.S,1)
    try
        Sm(j) = mean.S(j,find(~isnan(mean.S(j,:)),1,'first'));
        Tm(j) = mean.T(j,find(~isnan(mean.T(j,:)),1,'first'));
        Um(j) = mean.Vel(j,find(~isnan(mean.Vel(j,:)),1,'first'));
    catch
        Sm(j) = NaN;
        Tm(j) = NaN;
        Um(j) = NaN;
    end
end
[Tb,Sb,wb] = meltrate_gamm_dep_new(Tm',Sm',Um',topofile,0);

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

ax3 = subplot(3,1,3);
ax3_p = ax3.InnerPosition;
ax2_p = ax2.InnerPosition;
set(ax3,'InnerPosition',[ax3_p(1) ax3_p(2) ax2_p(3) ax3_p(4)]);
%wb(Ymesh(1,:)>=x_max) = 0;
plot(ax3,Ymesh,wb);
%set(ax3,'XTick',ax2,XTick);
xlabel(ax3,'Distance from GL (m)');
ylabel(ax3,'$\dot{m} (m/yr)$','interpreter','latex');
%xlim(ax3,[0 x_max]);
title(ax3,'Melt rate');

suptitle('Mean fields');
print(gcf,'-dpng','-r300',fullfile(pwd,'Mean_fields'));
%% Find the anamoly and co-variance field

for i=1:length(list_temp)
    
    T = squeeze(rdmds(list_temp(i).name(1:end-5)));
    S = squeeze(rdmds(list_salt(i).name(1:end-5)));
    S(S==0) = NaN;
    T(isnan(S)) = NaN;
    
    T(:,1:Z_index) = [];
    S(:,1:Z_index) = [];
    
    T = T(1:ver_skip:end,1:hor_skip:end);
    S = S(1:ver_skip:end,1:hor_skip:end);
    
    T_anom = T - mean.T;
    %T_anom_1 = T_anom'; 
    S_anom = S - mean.S;
    %S_anom_1 = S_anom';
    
    % Finding the unbiased co-variance field
    cov_T = (T_anom(:) * T_anom(:)')/(prod(size(T_anom)) - 1);
    cov_S = (S_anom(:) * S_anom(:)')/(prod(size(S_anom)) - 1);
    
    var_T = reshape(diag(cov_T),size(T_anom));
    var_S = reshape(diag(cov_S),size(S_anom));
    
end