%% This code will plot the temperature, salinity, freshwater flux,
% heat flux, melt rate and the associated density field for all the time
% steps corresponding to the 2D fields. Linear interpolation is used to
% interpolate the interface properties.
% Comment the next three lines if using with the shell script
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%% Output file directory and output files
file_dir = pwd;
list_T = dir(fullfile(file_dir,'T*.data'));
list_S = dir(fullfile(file_dir,'S.0*.data'));
list_surf = dir(fullfile(file_dir,'surf*.data'));
topofile = [file_dir filesep 'shelficeTopo.Lin.bin'];
delta_t = 4;
day2s = 24 * 3600;
day2hr = 24;
a = -5.73 * 10^-2;
b = 9.39 * 10^-2;
c = -7.53 * 10^-8;
g = 9.81;
rhoref = 1025;
X_max = 4000;
nx = max(size(squeeze(rdmds('DXC'))));
Y = cumsum(squeeze(rdmds('DYC')));
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);

%% Pressure estimation (Hydrostatic)
H = readbin(topofile,[1 nx])';
Pb = abs(H) * rhoref * g;

%% Timestamps of 2D fields
for i=1:length(list_S)
    time_2D(i) = str2num(list_S(i).name(end-10:end-5)) * delta_t;
end

%% Timestamps of the interface properties
for i=1:length(list_surf)
    time_interface(i) = str2num(list_surf(i).name(end-10:end-5)) * delta_t;
end

%% Extracting the interface properties
for i=1:length(list_surf)
    Data_2D = squeeze(rdmds(list_surf(i).name(1:end-5)));
    Um(i,:) = Data_2D(:,end);
    heat_flux(i,:) = Data_2D(:,9);
    fw_flux(i,:) = Data_2D(:,8);
end

%% Interpolating in time
for i=1:size(Um,2)
    Um_i(:,i) = interp1(time_interface,Um(:,i),time_2D);
    heat_flux_i(:,i) = interp1(time_interface,heat_flux(:,i),time_2D);
    fw_flux_i(:,i) = interp1(time_interface,fw_flux(:,i),time_2D);
end

time_2D = time_2D/day2s;
figure('units','normalized','outerposition',[0 0 1 1]);
%% Calculating the temperature, salinity and melt rate at the ice-shelf base
for i=1:length(list_S)
    clf;
    S = squeeze(rdmds(list_S(i).name(1:end-5)));
    T = squeeze(rdmds(list_T(i).name(1:end-5)));
    S(S==0) = NaN;
    S(Ymesh'>X_max) = NaN;
    T(isnan(S)) = NaN;
    
    for j=1:size(S,1)
        try
            Sm(i,j) = S(j,find(~isnan(S(j,:)),1,'first'));
            Tm(i,j) = T(j,find(~isnan(T(j,:)),1,'first'));
            %Um(i,j) = Vel(j,find(~isnan(Vel(j,:)),1,'first'));
        catch
            Sm(i,j) = NaN;
            Tm(i,j) = NaN;
            %Um(i,j) = NaN;
        end
    end
    
    [Tb(i,:),Sb(i,:),wb(i,:)] = meltrate_gamm_dep_new(Tm(i,:)',Sm(i,:)',Um(i,:)',topofile,0);
    
    Tf = a * Sb(i,:) + b + c * Pb';
    ax1 = subplot(2,3,1);
    contourf(ax1,Ymesh',Zmesh',S,'LineColor','none');colorbar;colormap(ax1,'jet');caxis(ax1,[34 35]);
    xlabel(ax1,'Distance from GL (m)');
    ylabel(ax1,'Depth (m)');
    title(ax1,'Salinity');
    
    ax2 = subplot(2,3,2);
    contourf(ax2,Ymesh',Zmesh',T,'LineColor','none');colorbar;colormap(ax2,'jet');caxis(ax2,[1.5 2]);
    xlabel(ax2,'Distance from GL (m)');
    ylabel(ax2,'Depth (m)');
    title(ax2,'Temperature');
    
    ax3 = subplot(2,3,3);
    plot(ax3,Ymesh,heat_flux_i(i,:));
    xlabel(ax3,'Distance from GL (m)');
    ylabel(ax3,'heat flux (W/m^2)');
    title(ax3,'Heat flux');
    
    ax4 = subplot(2,3,4);
    plot(ax4,Ymesh,fw_flux_i(i,:));
    xlabel(ax4,'Distance from GL (m)');
    ylabel(ax4,'FW flux (kg/m^2/s)');
    title(ax4,'FW flux');
    
    ax5 = subplot(2,3,5);
    plot(ax5,Ymesh,wb(i,:));
    xlabel(ax5,'Distance from GL (m)');
    ylabel(ax5,'$\dot{m} (m/yr)$','interpreter','latex');
    title(ax5,'Melt');
    
    ax6 = subplot(2,3,6);
    plot(ax6,Ymesh,Tb(i,:));
    xlabel(ax6,'Distance from GL (m)');
    ylabel(ax6,'Temp. at Base (C)');
    
    day = floor(time_2D(i));
    hours = (time_2D(i) - day) * day2hr;
    if(hours>0)
        time_vec = sprintf('%d day %d hours',day,hours);
    else
        time_vec = sprintf('%d day',day);
    end
    suptitle(sprintf('Time = %s',time_vec));
    print(gcf,'-dpng','-r300',fullfile(pwd,sprintf('Time_stamp_%04d',i)));
end
    


