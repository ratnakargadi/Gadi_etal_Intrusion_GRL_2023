%% This code will plot the heat flux time series at the three different 
% locations.
% COMMENT THE NEXT THREE LINES IF CALLING THROUGH BASH SCRIPT
%clear all;
%clc;
%close all;

%% 
addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
GZ = 6;
km2m = 10^3;
list_T = dir(fullfile(pwd,'T.*data'));
list_S = dir(fullfile(pwd,'S.*data'));
YC = cumsum(squeeze(rdmds('DYC')))/km2m - GZ;
dz = squeeze(rdmds('DRF'));
loc_want = [0];
delta_t = 0.3;
day2s = 3600 * 24;
T_freeze = -2.2;
nskip = 2;
%t_start = 12;
%t_end = 24;
n = t_end/12;

ind = find(YC<-0.1);

count = 1;
for i=1:nskip:length(list_T)
    code_T(count) = str2num(list_T(i).name(end-11:end-5)) * delta_t/day2s;
    S = squeeze(rdmds(list_S(i).name(1:end-5)));
    T = squeeze(rdmds(list_T(i).name(1:end-5)));
    %V = squeeze(rdmds(list_V(i).name(1:end-5)));
    
    T(S==0) = NaN;
    for j=1:length(ind)
        loc_val = find(~isnan(T(ind(j),:)));
        if(~isempty(loc_val))
            Tf(count,j) = T(ind(j),min(loc_val)) - T_freeze;
        else
            Tf(count,j) = -T_freeze;
        end
    end
    count = count + 1; 
   
end

time_ind = find(code_T*24>=t_start&code_T*24<t_end);

C = {'r','g','b','k','m','c'};
count = 1;

figure('units','normalized','outerposition',[0 0 1 1]);
for j=1:5:length(time_ind)
    plot(YC(ind)+GZ,Tf(time_ind(j),:),'Color',C{count},'LineWidth',2);
    hold on;
    L{count} = num2str(code_T(time_ind(j))*24,'%2.2f');
    count = count + 1;
end
LL = legend(L);
title (LL,'Time (hrs)');
LL.Title.Visible = 'on';
LL.FontSize = 16;
LL.FontWeight = 'bold';
LL.Location = 'Bestoutside';
xlabel('Distance from LF boundary (km)');
ylabel('T - T_f (\circC)');
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,sprintf('TF_%s',num2str(n,'%02d'))));

