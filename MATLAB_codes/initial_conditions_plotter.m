%% This code will plot the initial conditions for the simulation
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
nlevels = 200;
%% Mesh data
Y = cumsum(squeeze(rdmds('DYC')));
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);

S = squeeze(rdmds('S',0));
S(S==0) = NaN;

T = squeeze(rdmds('T',0));
T(isnan(S)) = NaN;

figure('units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(2,1,1);
contourf(ax1,Ymesh',Zmesh',S,nlevels,'LineStyle','none');
colormap(ax1,'jet');
colorbar;
%grid(ax1,'on');
xlabel(ax1,'Distance from GL (km)','FontSize',16,'FontWeight','bold');
ylabel(ax1,'Depth (m)','FontSize',16,'FontWeight','bold');
title(ax1,'Salinity','FontSize',16,'FontWeight','bold');

ax2 = subplot(2,1,2);
contourf(ax2,Ymesh',Zmesh',T,nlevels,'LineStyle','none');
colorbar;
xlabel(ax2,'Distance from GL (km)','FontSize',16,'FontWeight','bold');
ylabel(ax2,'Depth (m)','FontSize',16,'FontWeight','bold');
title(ax2,'Temperature','FontSize',16,'FontWeight','bold');

sgtitle('Initial Conditions','FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'IC'));