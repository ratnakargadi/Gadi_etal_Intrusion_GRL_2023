%% This code will plot the salinity, temperature and melt rate for three
% different time stamps
% Comment the next three lines after debugging if using in shell script
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
addpath(genpath('/nfspool-0/home/rgadi/MITgcm/gsw_matlab_v3_06_15'));
%% File directory and time stamp details
file_dir = pwd;
ver = 2;
tstart = 0;
delta_t = 2;
nlevel = 200;
g = 9.8;
rho_w = 1023;
%twant = [0.25;0.5;0.75]; % In days from the beginning data
%twant = [0;0.5;0.9];
twant = [0;2.5;5] + 10;
code_output = 10800;
day2s = 24 * 3600;
day2hr = 24;
Sfill = NaN;
Tfill = NaN;
x_max = 11000;
melt_x = 600;
Pa_dbar = 10^4;
km2m = 1000;
% Sfill = NaN;
% Tfill = NaN;
ymin_lim = -500;
ymax_lim = -300;
listsal = dir(fullfile(file_dir,'S.0*.data'));
listtemp = dir(fullfile(file_dir,'T*.data'));

for i=1:length(listsal)
    codet(i) = str2num(listsal(i).name(end-10:end-5))/day2s * delta_t;
end

% topofile = [file_dir filesep 'shelficeTopo.Lin.bin'];
% nx = length(squeeze(rdmds('DXC')));
% H = readbin(topofile,[1 nx]);
% P = rho_w * g * H/Pa_dbar;

Y = cumsum(squeeze(rdmds('DYC')));
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);
P = rho_w * g * abs(Zmesh)'/Pa_dbar;
rho_min = 10^5;rho_max = -10^5;

%% Using seawater toolbox to calculate density
for i=1:length(listsal)
    S = squeeze(rdmds(listsal(i).name(1:end-5)));
    T = squeeze(rdmds(listtemp(i).name(1:end-5)));
    S(S==0) = NaN;
    T(isnan(S)) = NaN;
    
    rho = gsw_rho(S,T,P);
    rho(end,:) = NaN;
    rho_min = min(rho_min,min(rho(:)));
    rho_max = max(rho_max,max(rho(:)));
    
    rho_3D(:,:,i) = rho;
    disp(sprintf('%d out of %d timesteps',i,length(listsal)));
end

%% Plotting the seawater density
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(listsal)
    clf;
    contourf(Ymesh'/km2m,Zmesh',squeeze(rho_3D(:,:,i)),nlevel,'LineStyle','none');
    colorbar;
    colormap('jet');
    xlabel('Distance from GL (km)');
    ylabel('Depth (m)');
    caxis([rho_min rho_max]);
    ylim([ymin_lim ymax_lim]);
    day = floor(codet(i));
    hour = codet(i) - day;
    if(hour>0)
        time_vec = sprintf('%d Day %d hours',day,hour*day2hr);
    else
        time_vec = sprintf('%d Day',day);
    end
    suptitle(sprintf('Density at Time = %s',time_vec));
    print(gcf,'-dpng','-r300',fullfile(file_dir,sprintf('dens_field_%04d',i)));
end


    