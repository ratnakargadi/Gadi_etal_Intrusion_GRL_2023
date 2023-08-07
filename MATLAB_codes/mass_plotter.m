clear all;
clc;
close all;

%%
addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
D = squeeze(rdmds('surfDiag',NaN));
Y = cumsum(squeeze(rdmds('DYC')));
delta_t = 1;
hours2s = 3600;
km2m = 10^3;
rho_w = 1028;
amp = 1;

shi_mass = squeeze(D(:,12,:));
t_start = 0;
t_end = 12;
list = dir(fullfile(pwd,'surfDiag.*data'));

for i=1:length(list)
    code_T(i) = str2num(list(i).name(end-11:end-5)) * delta_t/hours2s;
end

t_want = [3 6];

for j=1:length(t_want)
    ind(j) = min(find(abs(t_want(j) - code_T)==min(abs(t_want(j) - code_T))));
end

shi_mass(shi_mass==0) = NaN;
f = abs(shi_mass(:,ind(2)) - shi_mass(:,ind(1)))/rho_w;
ind_nan = find(isnan(f));
f(ind_nan(2)-1) = NaN;

Y = Y/km2m - 6;
loc_want = [-5 5 10];

for j=1:length(loc_want)
    ind(j) = find(abs(loc_want(j)-Y)==min(abs(loc_want(j)-Y)));
end
c = [0 0 0;1 0 1;0 1 0];
sz = 50; 

ind_nan = find(~isnan(shi_mass(:,ind(1))));
base_shape_1 = shi_mass(ind_nan,ind(1))/rho_w;
YY_1 = [Y(ind_nan) flip(Y(ind_nan))];
clear ind_nan;
ind_nan = find(~isnan(shi_mass(:,ind(2))));
base_shape_2 = shi_mass(ind_nan,ind(2))/rho_w;
YY_2 = [Y(ind_nan) flip(Y(ind_nan))];
clear ind_nan;
surf_elev_1 = base_shape_1/9.8;
surf_elev_2 = base_shape_2/9.8;

shape_1 = [-base_shape_1;flip(surf_elev_1)];
shape_2 = [-base_shape_2;flip(surf_elev_2)];

colororder({'k','r'})
plot(YY_1,shape_1,'b','LineWidth',2)
hold on;
plot(YY_2,shape_2,'k','LineWidth',2);
hold on;
plot([-6 -6],[-500 100],'--r','LineWidth',2);
hold on;
plot([0 0],[-500 100],'--b','LineWidth',2);
%hold on;
%plot([10 10],[-500 0],'--g','LineWidth',2);
hold on;
scatter(loc_want,[-shi_mass(ind(:),1)]/rho_w,sz,c,"filled");
hold off;
xticks([-6 0 10 20 30 40 50]);
xlim([-7 60]);
xlabel('Distance from GL (km)');
ylabel('Depth (m)');
yyaxis right
plot(Y,f,'r','LineWidth',2);
ylabel('Ice tidal motion (m)');
set(gca,'FontSize',18,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'Figure_2a'));


figure;
ind_t = find(code_T>=t_start&code_T<=t_end);
plot(code_T,-(shi_mass(ind(1),:)-mean(shi_mass(ind(1),ind_t)))/rho_w,'k','LineWidth',2);
hold on;
plot(code_T,-(shi_mass(ind(2),:)-mean(shi_mass(ind(2),ind_t)))/rho_w,'m','LineWidth',2);
hold on;
plot(code_T,-(shi_mass(ind(3),:)-mean(shi_mass(ind(3),ind_t)))/rho_w,'g','LineWidth',2);
ylabel('Ice tidal motion (m)');
xlabel('Time (hrs)');
xticks([0 3 6 9 12]);
xlim([t_start t_end]);
set(gca,'FontSize',18,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'Figure_2b'));


