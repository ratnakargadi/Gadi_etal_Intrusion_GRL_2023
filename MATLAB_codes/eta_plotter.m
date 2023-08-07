%% This code will plot the time series of eta. Useful to establish a remeshing
% frequency for shelf 2D remesh code
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
list_eta = dir(fullfile(pwd,'Eta*.data'));
list_FW = dir(fullfile(pwd,'SHICE_f*.data'));
list_heat = dir(fullfile(pwd,'SHICE_h*.data'));

delta_t = 25;
Y = squeeze(rdmds('YC'));
if(length(unique(Y))==1)
    Y = squeeze(rdmds('XC'));
end

eta_min = 10^3;
eta_max = -10^3;
FW_min = 10^9;
FW_max = -10^9;
heat_min = 10^9;
heat_max = -10^9;

for i=1:length(list_eta)
    eta = squeeze(rdmds(list_eta(i).name(1:end-5)));
    eta_min = min(eta_min,min(eta));
    eta_max = max(eta_max,max(eta));
    try
        FW = squeeze(rdmds(list_FW(i).name(1:end-5)));
        heat = squeeze(rdmds(list_heat(i).name(1:end-5)));
    
        FW_min = min(FW_min,min(FW));
        FW_max = max(FW_max,max(FW));
        heat_min = min(heat_min,min(heat));
        heat_max = max(heat_max,max(heat));
    end
end

figure('units','normalized','outerposition',[0 0 1 1]);

for i=1:length(list_eta)
    try
        clf;
        eta = squeeze(rdmds(list_eta(i).name(1:end-5)));
        FW = squeeze(rdmds(list_FW(i).name(1:end-5)));
        heat = squeeze(rdmds(list_heat(i).name(1:end-5)));
    
        ax1 = subplot(1,3,1);
        plot(ax1,Y,eta);
        xlabel(ax1,'Longitude (\circ)');
        ylabel(ax1,'\eta');
        ylim(ax1,[eta_min eta_max]);
    
        ax2 = subplot(1,3,2);
        plot(ax2,Y,FW);
        xlabel(ax2,'Longitude (\circ)');
        ylabel(ax2,'SHICE_fwFlux','Interpreter','none');
        ylim(ax2,[FW_min FW_max]);
    
        ax3 = subplot(1,3,3);
        plot(ax3,Y,heat);
        xlabel(ax3,'Longitude (\circ)');
        ylabel(ax3,'SHICE_heatFlux','Interpreter','none');
        ylim(ax3,[heat_min heat_max]);
    
        suptitle(sprintf('Time = %d s',delta_t * str2num(list_eta(i).name(end-7:end-5))));
    end
    
    print(gcf,'-dpng','-r300',fullfile(pwd,sprintf('ETA_%s',num2str(i,'%04d'))));
end