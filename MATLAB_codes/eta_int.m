%% This code plots the instantaneous velocity along the GL
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%%
list_eta = dir(fullfile(pwd,'surfDiag*.data'));
deg2km = 111.4;
delta_t = 1;
s2hr = 3600;
dist_from_ice_f = 25;

for i=1:length(list_eta)
    se = squeeze(rdmds(list_eta(i).name(1:end-5)));
    time(i) = str2num(list_eta(i).name(end-14:end-5)) * delta_t/s2hr;
    GL_j = squeeze(se(:,7,:));
    GL(i) = GL_j(1,1);
    eta(i,:) = squeeze(se(:,1,:));
end
hFac = squeeze(rdmds('hFacC'));
YC = squeeze(rdmds('YC'));
hIce = readbin('shelficeTopo.Lin.bin',[1 length(YC)]);
ind_f = max(find(hIce==0,2));
ind_f1 = max(find(YC<=YC(ind_f)-dist_from_ice_f/deg2km));
ind_f2 = max(find(YC<=YC(ind_f)-dist_from_ice_f*2/deg2km));

%% Integrating eta
for i=1:length(list_eta)
    eta_int_calc(i) = trapz(YC(ind_f:end)*deg2km,eta(i,ind_f:end))/trapz(YC(ind_f:end)*deg2km);
    eta_int_calc_n(i) = sum(eta(i,ind_f:end))/(length(YC) - ind_f);
    eta_front(i) = eta(i,ind_f);
    eta_last_ice_shelf(i) = eta(i,ind_f-1);
    eta_sp(i) = eta(i,ind_f1);
    eta_spp(i) = eta(i,ind_f2);
end

%% Plotting the integrated free surface elevation
plot(time,eta_int_calc_n,'k','LineWidth',2);
xlabel('Time (hr)');
ylabel({'Integrated free surface','elevation (\eta)'});
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'Integrated_eta'));

figure('Units','normalized','OuterPosition',[0 0 1 1]);
subplot(4,1,1);
plot(time,eta_front,'k','LineWidth',2);
xlabel('Time (hr)');
ylabel({'Free surface','1st cell ','after ice-front'});
set(gca,'FontSize',16,'FontWeight','bold');

subplot(4,1,2);
plot(time,eta_last_ice_shelf,'k','LineWidth',2);
xlabel('Time (hr)');
ylabel({'Free surface','at ice-front'});
set(gca,'FontSize',16,'FontWeight','bold');

subplot(4,1,3);
plot(time,eta_sp,'k','LineWidth',2);
xlabel('Time (hr)');
str = [sprintf('%d km',dist_from_ice_f) ' from'];
ylabel({'Free surface at',str,' ice-front'});
set(gca,'FontSize',16,'FontWeight','bold');

subplot(4,1,4);
plot(time,eta_spp,'k','LineWidth',2);
xlabel('Time (hr)');
str = [sprintf('%d km',2*dist_from_ice_f) ' from'];
ylabel({'Free surface at',str,' ice-front'});
set(gca,'FontSize',16,'FontWeight','bold');

print(gcf,'-dpng','-r300',fullfile(pwd,'eta_three_points'));
