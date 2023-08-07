%% This code will use the freshwater flux data and divide it with density of
% freshwater to get the melt rate. This will be converted from m/s to m/yr
% and plotted
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%% File directory with the list of all the melt rate files
file_dir = pwd;
list = dir(fullfile(file_dir,'SHICE_fwFlux*.data'));


%% Parameters
km2m = 10^3;
delta_t = 1;
GZ = 6;
day2s = 24 * 3600;
day2hr = 24;
rho_w = 10^3;
rho_I = 917;
y2s = 365 * 24 * 3600;
km2m = 10^3;
ncycle_want = 7;
tidal_period = 12;
time_min = 0;
time_max = tidal_period*ncycle_want;
tvec = time_min:tidal_period:time_max;

%% Mesh data
Y = cumsum(squeeze(rdmds('DYC')));
dY = squeeze(rdmds('YC'));
Y = Y - GZ*km2m;

%disp(length(list));
%error('A');
%figure;
inds = find(Y/10^3<-0.1);

for j=1:length(tvec)-1
    Fw_fl = zeros(1,length(Y));
    count = 0;
        for i=1:length(list)
%    i=length(list);
        code_T = str2num(list(i).name(end-11:end-5)) * delta_t/day2s;
        Fw_flux(i,:) = y2s * squeeze(rdmds(list(i).name(1:end-5)))/rho_w;
        if(24*code_T>=tvec(j)&&24*code_T<=tvec(j+1))
            Fw_fl = Fw_fl + Fw_flux(i,:);
            count = count + 1;
        end
    
        melt_avg(i) = sum(Fw_flux(i,:) .* dY)/(Y(end)); 
        end

    melt_tidal_avg = abs(Fw_fl)/count;
    melt_int(j) = trapz(Y(inds)/10^3-0.2,melt_tidal_avg(inds)) * rho_I/rho_w/(GZ);
 
end

plot(1:length(tvec)-1,melt_int,'k','LineWidth',2);
ylabel('Mean melt rate (m/yr)');
xlabel('Tidal cycle');
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'mean_melt_cycle'));

save(sprintf('Cycle_data_%d.mat',num2str(GZ,'%02d')),'tvec','melt_int','GZ');
