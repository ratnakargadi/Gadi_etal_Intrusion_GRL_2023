%% This code computes the total oceanic heat content to establish the 
% convergence to steady state
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%% Parameters
file_dir = pwd;
list_T = dir(fullfile(file_dir,'T.0*data'));
list_S = dir(fullfile(file_dir,'S.0*data'));
list_V = dir(fullfile(file_dir,'V.0*data'));
list_W = dir(fullfile(file_dir,'W.0*data'));
list_fwflx = dir(fullfile(file_dir,'SHICE_fwFlux*.data')); 
list_heat = dir(fullfile(file_dir,'SHICE_h*.data'));

dY = (squeeze(rdmds('DYC')));
dZ = -squeeze(rdmds('DRF'));
Z = cumsum(dZ);
Y = cumsum(dY);

delta_t = 1.5;
day2s = 24 * 3600;
day2hr = 24;
y2s = 365 * 24 * 3600;

[Ymesh,Zmesh] = meshgrid(Y,Z);
Area = dY' * abs(dZ)';
rho_w = 1023; % later use the seawater toolbox to compute the density. Use reference pressure as
% the pressure at the base of the ice-shelf
Cp = 3974;
%% Extracting the temperature field
for i=1:length(list_T)
    S = squeeze(rdmds(list_S(i).name(1:end-5)));
    S(S==0) = NaN;
    T = squeeze(rdmds(list_T(i).name(1:end-5)));
    T(isnan(S)) = 0;
    V = squeeze(rdmds(list_V(i).name(1:end-5)));
    V(isnan(S)) = 0;
    W = squeeze(rdmds(list_W(i).name(1:end-5)));
    
    try
        Fw_flux = y2s * squeeze(rdmds(list_fwflx(i).name(1:end-5)))/rho_w;
        heat_f = squeeze(rdmds(list_heat(i).name(1:end-5)));
        Fw_flux(find(Fw_flux<=0)) = -Fw_flux(find(Fw_flux<=0));
        melt_avg(i) = sum(Fw_flux .* dY)/(Y(end));
        surf_heat_flux(i) = sum(heat_f .* dY);
        code_T_surf(i) = str2num(list_heat(i).name(end-11:end-5)) * delta_t/day2s;
    end
    Heat_content(i) = sum(T.*Area.*rho_w,'all') * Cp;
    KE(i) = sum(V.*V.*Area.*rho_w,'all')*1/2 + 1/2*sum(W.*W.*Area.*rho_w,'all');
    code_T(i) = str2num(list_S(i).name(end-11:end-5)) * delta_t/day2s;
    disp(i);
end

%% Plotting 
figure('units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(2,4,1);
plot(ax1,code_T,Heat_content,'LineWidth',2,'Color','k');
xlabel(ax1,'Time (days)','FontSize',16,'FontWeight','bold');
ylabel(ax1,'Heat Content (J)','FontSize',16,'FontWeight','bold');

ax2 = subplot(2,4,2);
plot(ax2,code_T,KE,'LineWidth',2,'Color','k');
xlabel(ax2,'Time (days)','FontSize',16,'FontWeight','bold');
ylabel(ax2,'Kinetic Energy (J)','FontSize',16,'FontWeight','bold');

ax3 = subplot(2,4,3);
plot(ax3,code_T_surf,melt_avg,'LineWidth',2,'Color','k');
xlabel(ax3,'Time (days)','FontSize',16,'FontWeight','bold');
ylabel(ax3,{'Spatially averaged';'Melt rate (m/yr)'},'FontSize',16,'FontWeight','bold');

ax4 = subplot(2,4,4);
plot(ax4,code_T_surf,surf_heat_flux,'LineWidth',2,'Color','k');
xlabel(ax4,'Time (days)','FontSize',16,'FontWeight','bold');
ylabel(ax4,{'Spatially averaged';'heat flux'},'FontSize',16,'FontWeight','bold');

t1 = (code_T(1:end-1) + code_T(2:end))/2;
ax5 = subplot(2,4,5);
plot(ax5,t1,diff(Heat_content),'LineWidth',2,'Color','k');
xlabel(ax5,'Time (days)','FontSize',16,'FontWeight','bold');
ylabel(ax5,'Diff. Heat Content (J)','FontSize',16,'FontWeight','bold');

ax6 = subplot(2,4,6);
plot(ax6,t1,diff(KE),'LineWidth',2,'Color','k');
xlabel(ax6,'Time (days)','FontSize',16,'FontWeight','bold');
ylabel(ax6,'diff. KE (J)','FontSize',16,'FontWeight','bold');

t2 = (code_T_surf(1:end-1) + code_T_surf(2:end))/2;
ax7 = subplot(2,4,7);
plot(ax7,t2,diff(melt_avg),'LineWidth',2,'Color','k');
xlabel(ax7,'Time (days)','FontSize',16,'FontWeight','bold');
ylabel(ax7,{'Diff. spatially averaged';'melt rate (m/yr)'},'FontSize',16,'FontWeight','bold');

ax8 = subplot(2,4,8);
plot(ax8,t2,diff(surf_heat_flux),'LineWidth',2,'Color','k');
xlabel(ax8,'Time (days)','FontSize',16,'FontWeight','bold');
ylabel(ax8,{'Diff. spatially averaged';'heat flux'},'FontSize',16,'FontWeight','bold');

sgtitle('Convergence');
print(gcf,'-dpng','-r300',fullfile(pwd,'Convergence'));
