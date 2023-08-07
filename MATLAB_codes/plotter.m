%% This code will plot the salinity, temperature and melt rate for three
% different time stamps
% Comment the next three lines after debugging if using in shell script
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%% File directory and time stamp details
file_dir = pwd;
ver = 3;
tstart = 0;
delta_t = 2;
nlevels = 2000;
%twant = [0.25;0.5;0.75]; % In days from the beginning data
%twant = [0;0.5;0.9];
twant = [0;2.5;5] + 15;
code_output = 10800;
day2s = 24 * 3600;
Sfill = NaN;
Tfill = NaN;
x_max = 11000;
melt_x = 600;
% Sfill = NaN;
% Tfill = NaN;
ymin_lim = -500;
ymax_lim = -300;
list = dir(fullfile(file_dir,'dynDiag*.data'));
listsal = dir(fullfile(file_dir,'S.0*.data'));
listtemp = dir(fullfile(file_dir,'T*.data'));
%listsurf = dir(fullfile(file_dir,'surfDiag*.data'));
listU = dir(fullfile(file_dir,'U*.data'));
listV = dir(fullfile(file_dir,'V*.data'));
listW = dir(fullfile(file_dir,'W*.data'));
listsurf = dir(fullfile(file_dir,'surfDiag*.data'));
for i=1:length(listsal)
    codet(i) = str2num(listsal(i).name(end-10:end-5));
end
for i=1:length(listsurf)
    codesurft(i) = str2num(listsurf(i).name(end-10:end-5));
end
deg2km = 111.4;
km2m = 1000;
topofile = 'shelficeTopo.Lin.bin';

%% Creating the iteration number/ time stamp number
timestamp = twant * day2s/delta_t + tstart;

for i=1:length(timestamp)
   ind = find(timestamp(i)==codet);
   inds = find(timestamp(i)==codesurft);
   if(isempty(ind))
       dist = abs(timestamp(i) - codet);
       ind = min(find(min(dist)==dist));
   end
   if(isempty(inds))
       dist = abs(timestamp(i) - codesurft);
       inds = min(find(min(dist)==dist));
   end
   index(i) = ind;
   indexs(i) = inds;
end

%Y = squeeze(rdmds('YC'));
Y = cumsum(squeeze(rdmds('DYC')));
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);

for i=1:length(index)
%     Data = squeeze(rdmds(list(index(i)).name(1:end-5)));
%     S = squeeze(Data(:,:,2+3));
%     T = squeeze(Data(:,:,2+2));
    S = squeeze(rdmds(listsal(index(i)).name(1:end-5)));
    T = squeeze(rdmds(listtemp(index(i)).name(1:end-5)));
    U = squeeze(rdmds(listU(index(i)).name(1:end-5)));
    V = squeeze(rdmds(listV(index(i)).name(1:end-5)));
    W = squeeze(rdmds(listW(index(i)).name(1:end-5)));
    %Vel = sqrt(U.^2 + V.^2 + W.^2);
    Vel = sqrt(U.^2 + V.^2);
    Vel(Vel==0) = NaN;
    S(S==0) = NaN;
    T(T==0) = NaN;
    
    clear Data;
    Data_2D = squeeze(rdmds('surfDiag',codesurft(indexs(i))));
    %Data_2D = squeeze(rdmds('surfDiag',433350));
    %Uc = squeeze(Data_2D(:,end));
    %% LOOK AT THIS SECTION LATER
%     try
%         Um(i,:) = squeeze(Data_2D(:,end));
%     catch
%         Um(i,:) = Um(i-1,:);
%     end
%     
    for j=1:size(S,1)
        try
            Sm(i,j) = S(j,find(~isnan(S(j,:)),1,'first'));
            Tm(i,j) = T(j,find(~isnan(T(j,:)),1,'first'));
            Umm(i,j) = Vel(j,find(~isnan(Vel(j,:)),1,'first'));
        catch
            Sm(i,j) = NaN;
            Tm(i,j) = NaN;
            Umm(i,j) = NaN;
        end
    end
    Um(i,:) = squeeze(Data_2D(:,end));
    S(isnan(S)) = Sfill;
    T(isnan(T)) = Tfill;
    S(end,:) = NaN;
    T(end,:) = NaN;
    Salt(:,:,i) = S;
    Temp(:,:,i) = T;
   [Tb(i,:),Sb(i,:),wb(i,:)] = meltrate_gamm_dep_new(Tm(i,:)',Sm(i,:)',Um(i,:)',topofile,0);
   %[Tb(i,:),Sb(i,:),wb(i,:)] = meltrate_gamm_dep_new(Tm(i,:)',Sm(i,:)',Umm(i,:)',topofile,1);
end

figure('units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(3,3,1);
contourf(ax1,Ymesh,Zmesh,squeeze(Salt(:,:,1))',nlevels,'LineColor','none');colorbar;colormap(ax1,'jet');
Salt_min = min(Salt(:)); Salt_max = max(Salt(:));
%Salt_min = 33; Salt_max = 35.5;
ylim([ymin_lim ymax_lim]);
caxis(ax1,[Salt_min Salt_max]);
xlabel(ax1,'Distance from GL (m)');
ylabel(ax1,'Depth (m)');
title(ax1,sprintf('Salinity after %d days',twant(1)));

ax2 = subplot(3,3,4);
contourf(ax2,Ymesh,Zmesh,squeeze(Temp(:,:,1))',nlevels,'LineColor','none');colorbar;colormap(ax2,'jet');
Temp_min = min(Temp(:));Temp_max = max(Temp(:));
ylim([ymin_lim ymax_lim]);
caxis(ax2,[Temp_min Temp_max]);
xlabel(ax2,'Distance from GL (m)');
ylabel(ax2,'Depth (m)');
title(ax2,sprintf('Temp. after %d days',twant(1)));

ax3 = subplot(3,3,7);
plot(ax3,Ymesh,wb(1,:));
xlabel(ax3,'Distance from GL (m)');
ylabel(ax3,'$\dot{m} (m/yr)$','interpreter','latex');
title(ax3,sprintf('Melt rate after %d days',twant(1)));
xlim(ax3,[0 x_max]);

h = axes('Position',[ax3.Position(1)+0.04 ax3.Position(2)+0.1 0.05 0.1]);
plot(h,Ymesh,wb(1,:));
xlabel(h,'Distance from GL (m)');
ylabel(h,'$\dot{m} (m/yr)$','interpreter','latex');
xlim(h,[0 melt_x]);

ax4 = subplot(3,3,2);
contourf(ax4,Ymesh,Zmesh,squeeze(Salt(:,:,2))',nlevels,'LineColor','none');colorbar;colormap(ax4,'jet');
caxis(ax4,[Salt_min Salt_max]);
ylim([ymin_lim ymax_lim]);
xlabel(ax4,'Distance from GL (m)');
ylabel(ax4,'Depth (m)');
title(ax4,sprintf('Salinity after %d days',twant(2)));

ax5 = subplot(3,3,5);
contourf(ax5,Ymesh,Zmesh,squeeze(Temp(:,:,2))',nlevels,'LineColor','none');colorbar;colormap(ax5,'jet');
caxis(ax5,[Temp_min Temp_max]);
ylim([ymin_lim ymax_lim]);
xlabel(ax5,'Distance from GL (m)');
ylabel(ax5,'Depth (m)');
title(ax5,sprintf('Temp. after %d days',twant(2)));

ax6 = subplot(3,3,8);
plot(ax6,Ymesh,wb(2,:));
xlabel(ax6,'Distance from GL (m)');
ylabel(ax6,'$\dot{m} (m/yr)$','interpreter','latex');
title(ax6,sprintf('Melt rate after % days',twant(2)));
xlim(ax6,[0 x_max]);

h = axes('Position',[ax6.Position(1)+0.04 ax6.Position(2)+0.1 0.05 0.1]);
plot(h,Ymesh,wb(2,:));
xlabel(h,'Distance from GL (m)');
ylabel(h,'$\dot{m} (m/yr)$','interpreter','latex');
xlim(h,[0 melt_x]);

ax7 = subplot(3,3,3);
contourf(ax7,Ymesh,Zmesh,squeeze(Salt(:,:,3))',nlevels,'LineColor','none');colorbar;colormap(ax7,'jet');
caxis([Salt_min Salt_max]);
ylim([ymin_lim ymax_lim]);
xlabel(ax7,'Distance from GL (m)');
ylabel(ax7,'Depth (m)');
title(ax7,sprintf('Salinity after %d days',twant(3)));

ax8 = subplot(3,3,6);
contourf(ax8,Ymesh,Zmesh,squeeze(Temp(:,:,3))',nlevels,'LineColor','none');colorbar;colormap(ax8,'jet');
ylim([ymin_lim ymax_lim]);
caxis(ax2,[Temp_min Temp_max]);
xlabel(ax8,'Distance from GL (m)');
ylabel(ax8,'Depth (m)');
title(ax8,sprintf('Temp. after %d days',twant(3)));

ax9 = subplot(3,3,9);
plot(ax9,Ymesh,wb(3,:));
xlabel(ax9,'Distance from GL (m)');
ylabel(ax9,'$\dot{m} (m/yr)$','interpreter','latex');
title(ax9,sprintf('Melt rate after %d days',twant(3)));
xlim(ax9,[0 x_max]);

h = axes('Position',[ax9.Position(1)+0.04 ax9.Position(2)+0.1 0.05 0.1]);
plot(h,Ymesh,wb(3,:));
xlabel(h,'Distance from GL (m)');
ylabel(h,'$\dot{m} (m/yr)$','interpreter','latex');
xlim(h,[0 melt_x]);

print(gcf,'-dpng','-r300',fullfile(pwd,sprintf('inst_hr_new_%d',ver)));