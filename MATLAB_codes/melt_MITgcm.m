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
%plot_dir = [pwd filesep 'MELT'];
%work_dir = '/DFS-L/DATA/rignot/rgadi/08302022/run_16'; 
work_dir = pwd;
plot_dir = [work_dir filesep 'MELT'];
if(~exist(plot_dir))
  mkdir(plot_dir);
end

%% Parameters
delta_t = 1;
day2s = 24 * 3600;
day2hr = 24;
rho_w = 10^3;
y2s = 365 * 24 * 3600;
km2m = 10^3;
nskip = 1;
%% Mesh data
Y = cumsum(squeeze(rdmds('DYC')));
GL_ini = Y(709);
Y = Y - 5000;
dY = squeeze(rdmds('DYC'));
%disp(length(list));
%error('A');
figure;
for i=1:nskip:length(list)
%    i=length(list);
    code_T = str2num(list(i).name(end-11:end-5)) * delta_t/day2s;
    Fw_flux(i,:) = y2s * squeeze(rdmds(list(i).name(1:end-5)))/rho_w;
    Fw_flux(:,find(Fw_flux(i,:)<=0)) = -Fw_flux(:,find(Fw_flux(i,:)<=0));
    
    day = floor(code_T);
    hour = code_T - day;
    if(hour>0)
        time_vec = sprintf('%d Day %d hours',day,hour*day2hr);
    else
        time_vec = sprintf('%d Day',day);
    end
    melt_avg(i) = sum(Fw_flux(i,:) .* dY)/(Y(end)); 
    clf;
    plot(Y/km2m,Fw_flux(i,:),'LineWidth',2);
    xlabel('Distance from GL (km)','FontSize',12,'FontWeight','bold');
    ylabel('$\dot{m} (m/yr)$','interpreter','latex','FontSize',12,'FontWeight','bold');
    %ylim([0 8]);
    title(sprintf('Time = %s',time_vec),'FontSize',12,'FontWeight','bold');
    %ylim([0 30]);
   % p = get(gca, 'Position');
   % h = axes('Parent', gcf, 'Position', [p(1)+0.05 p(2)+0.4 p(3)-0.55 p(4)-0.45]);
   % plot(h,Y/km2m,Fw_flux(i,:),'LineWidth',2);
   % xlim(h,[0 2]);
    %ylim(h,[0 4]);
    xlim([-10 10]);
    pause(1);
    %print(gcf,'-dpng','-r300',fullfile(plot_dir,sprintf('Melt_rate_%04d',i)));
    tvec(i) = code_T;
end
    
save(fullfile(work_dir,'Melt_average.mat'),'melt_avg','tvec');
