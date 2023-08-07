clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));

Y = cumsum(squeeze(rdmds('DYC')));
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);

t0 = [0 80 560 960 1840 8800];
delta_t = 5;
ylim_1 = -225;ylim_2 = -220;

figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'Color','w');

for i=1:length(t0)
    S_2D = squeeze(rdmds('S',t0(i)));
    S_2D(S_2D==0) = NaN;
    
    ax = subplot(3,length(t0)/3,i);
    contourf(ax,Ymesh,Zmesh,S_2D');
    ylim(ax,[ylim_1 ylim_2]);
    colorbar;
    xlabel(ax,'Distance from GL (m)');
    ylabel(ax,'Depth (m)');
    caxis(ax,[33 35]);
    title(ax,sprintf('time = %d s',t0(i)*delta_t));
end