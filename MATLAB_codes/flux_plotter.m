%% This code will plot the heat flux time series at the three different 
% locations.
% COMMENT THE NEXT THREE LINES IF CALLING THROUGH BASH SCRIPT
clear all;
clc;
close all;

%% 
addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
GZ = 6;
km2m = 10^3;
list_T = dir(fullfile(pwd,'T.*data'));
list_S = dir(fullfile(pwd,'S.*data'));
list_V = dir(fullfile(pwd,'V.*data'));
list_dyn = dir(fullfile(pwd,'dynDiag.*data'));
list_fw = dir(fullfile(pwd,'SHICE_fw*.data'));
YC = cumsum(squeeze(rdmds('DYC')))/km2m - GZ;
dz = squeeze(rdmds('DRF'));
loc_want = [-GZ+1 -GZ/2 0];
delta_t = 0.3;
day2s = 3600 * 24;
T_freeze = -2;
rho_w = 1028;
y2s = 365 * 24 * 3600;

for j=1:length(loc_want)
    ind(j) = min(find(abs(YC-loc_want(j))==min(abs(YC-loc_want(j)))));
end

for i=1:length(list_T)
    code_T(i) = str2num(list_T(i).name(end-11:end-5)) * delta_t/day2s;
    S = squeeze(rdmds(list_S(i).name(1:end-5)));
    T = squeeze(rdmds(list_T(i).name(1:end-5)));
    V = squeeze(rdmds(list_V(i).name(1:end-5)));
    
    T(S==0) = NaN;
    for j=1:length(ind)
        Tf(i,j) = rho_w * sum((squeeze(T(ind(j),:))-T_freeze)'.*dz.*squeeze(V(ind(j),:)'),'omitnan');
        Mf(i,j) = rho_w * sum(dz.*squeeze(V(ind(j),:)'),'omitnan');
        %Melt(i,j) = fw(ind(j));
    end
end

for i=1:length(list_fw)
    fw = squeeze(rdmds(list_fw(i).name(1:end-5))) * y2s/rho_w;
    code_F(i) = str2num(list_fw(i).name(end-11:end-5)) * delta_t/day2s;
    for j=1:length(ind)
        Melt(i,j) = fw(ind(j));
    end
end

figure('units','normalized','outerposition',[0 0 1 1]);
plot(code_T*24,Tf(:,1),'b','LineWidth',2)
hold on
plot(code_T*24,Tf(:,2),'r','LineWidth',2)
hold on
plot(code_T*24,Tf(:,3),'k','LineWidth',2)
L = legend('GL - 5 km','GL - 3 km','GL');
L.FontSize = 16;
L.FontWeight = 'bold';
ylim([-1 1]*max(abs(Tf(:))*1.2));
L.Location = 'best';
xlabel('Time (hrs)');
ylabel('Heat flux (\rho*V*(T-T_f))');
%xlim([0 100]);
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'heat_flux_potential_1'));

figure('units','normalized','outerposition',[0 0 1 1]);
plot(code_T*24,Mf(:,1),'b','LineWidth',2)
hold on
plot(code_T*24,Mf(:,2),'r','LineWidth',2)
hold on
plot(code_T*24,Mf(:,3),'k','LineWidth',2)
L = legend('GL - 5 km','GL - 3 km','GL');
L.FontSize = 16;
L.FontWeight = 'bold';
ylim([-1 1]*max(abs(Mf(:)))*1.2);
L.Location = 'best';
xlabel('Time (hrs)');
ylabel('Mom (\rho*V)');
%xlim([0 100]);
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'Momentum_1'));

figure('units','normalized','outerposition',[0 0 1 1]);
plot(code_F*24,-Melt(:,1),'b','LineWidth',2);
hold on;
plot(code_F*24,-Melt(:,2),'r','LineWidth',2);
hold on;
plot(code_F*24,-Melt(:,3),'k','LineWidth',2);
%xlim([0 100]);
ylim([0 1]*max(abs(Melt(:)))*1.2);
L = legend('GL - 5 km','GL - 3km','GL');
L.Location = 'Best';
L.FontSize = 16;
L.FontWeight = 'bold';
xlabel('Time (hrs)');
ylabel('Melt rate (m/yr)');
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'Melt_time_plot_1'));
