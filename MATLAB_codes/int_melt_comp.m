%% This code will use the freshwater flux data and divide it with density of
% freshwater to get the melt rate. This will be converted from m/s to m/yr
% and plotted
%clear all;
%clc;
%close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%% File directory with the list of all the melt rate files
file_dir = pwd;
list = dir(fullfile(file_dir,'SHICE_fwFlux*.data'));
%plot_dir = [pwd filesep 'MELT'];
%work_dir = '/DFS-L/DATA/rignot/rgadi/08302022/run_16'; 
%work_dir = pwd;
%plot_dir = [work_dir filesep 'MELT'];
%if(~exist(plot_dir))
%  mkdir(plot_dir);
%end

%% Parameters
delta_t = 1;
day2s = 24 * 3600;
day2hr = 24;
rho_I = 917;
rho_w = 10^3;
y2s = 365 * 24 * 3600;
km2m = 10^3;
nskip = 1;
%time_min = 0;
%time_max = 12;
n = floor(time_max/12);
width = 20;
%% Mesh data
Y = cumsum(squeeze(rdmds('DYC')));
dY = squeeze(rdmds('YC'));
%disp(length(list));
%error('A');
%figure;
Fw_fl = zeros(1,length(Y));
count = 0;
for i=1:nskip:length(list)
%    i=length(list);
    code_T = str2num(list(i).name(end-11:end-5)) * delta_t/day2s;
    Fw_flux(i,:) = y2s * squeeze(rdmds(list(i).name(1:end-5)))/rho_w;
    if(24*code_T>=time_min&&24*code_T<=time_max)
        Fw_fl = Fw_fl + Fw_flux(i,:);
        count = count + 1;
    end
    %Fw_flux(:,find(Fw_flux(i,:)<=0)) = -Fw_flux(:,find(Fw_flux(i,:)<=0));
    
    day = floor(code_T);
    hour = code_T - day;
    if(hour>0)
        time_vec = sprintf('%d Day %d hours',day,hour*day2hr);
    else
        time_vec = sprintf('%d Day',day);
    end
    melt_avg(i) = sum(Fw_flux(i,:) .* dY)/(Y(end)); 
    tvec(i) = code_T;
end
    
%save(fullfile(work_dir,'Melt_average.mat'),'melt_avg','tvec');
melt_tidal_avg = abs(Fw_fl)/count;

%% Plot
%figure('units','normalized','outerposition',[0 0 1 1]);
%figure;
%GL_ini = Y(489);
GL_init = 6 * 10^3;
%y_max = 80;
Y = Y - GL_init;
%GL_ini = GL_ini/10^3;
%GL_ini = 0;
%plot(Y/10^3,melt_tidal_avg,'k','LineWidth',2);
%hold on;
%plot([GL_ini/10^3 GL_ini/10^3],[0 y_max],'--b','LineWidth',2);
%hold on;
%plot(-[GL_init/10^3 GL_init/10^3],[0 y_max],'--r','LineWidth',2);
%hold on;
%plot([GL_ini/10^3 GL_ini/10^3]-2.5*2*1,[0 y_max],'--g','LineWidth',2);
%xlabel('Distance from GL (km)');
%ylabel('Basal Melt Rate (m/yr)');
%T = title({'Basal Melt Rate averaged','on semi-diurnal tidal period'});
%T.FontSize = 16;
%T.FontWeight = 'bold';
%set(gca,'FontSize',16,'FontWeight','bold');
%xlim([-10 50]);
%print(gcf,'-dpng','-r300',fullfile(pwd,'melt_rate_tidally_averaged'));
%hold on;
%plot([4.9 4.9],[0 y_max],'c','LineStyle','--','LineWidth',2);
%hold on;
%plot([10 10],[0 y_max],'m','LineStyle','--','LineWidth',2);
%save(sprintf('tidal_average_%d_%d',GL_init/10^3,n),'Y','GL_init','melt_tidal_avg');
%load('Melt_converged.mat');
%hold on;
%plot(Y/10^3+4.9,melt_p,'c','LineWidth',2);
%L = legend('$\dot{m}$ (GL moving)','GL (low tide)','GL (high tide)','GL (2015)','Flexure zone','$\dot{m}$ (Fixed GL)','Interpreter','latex');
%L.FontSize = 16;
%L.FontWeight = 'bold';
%print(gcf,'-dpng','-r300',fullfile(pwd,'melt_rate_tidally_averaged'));

inds_GZ = find(Y/10^3<-0.1);
inds_FZ = find(Y/10^3>=-0.1&Y/10^3<=10);
inds_else = find(Y/10^3>10);

int_melt.GZ = trapz(Y(inds_GZ)/10^3,melt_tidal_avg(inds_GZ))*rho_w*width*km2m/10^9;
int_melt.FZ = trapz(Y(inds_FZ)/10^3,melt_tidal_avg(inds_FZ))*rho_w*width*km2m/10^9;
int_melt.else = trapz(Y(inds_else)/10^3,melt_tidal_avg(inds_else))*rho_w*width*km2m/10^9;

disp(sprintf('GZ melt = %4.2f Gt/yr',int_melt.GZ));
disp(sprintf('FZ melt = %4.2f Gt/yr',int_melt.FZ));
disp(sprintf('Melt elsewhere = %4.2f Gt/yr',int_melt.else));
