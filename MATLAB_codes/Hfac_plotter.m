%% This code will plot the 2D velocity fields time series under an ice-sheet
% Comment the next three lines after debugging
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%addpath(genpath('/DFS-L/DATA/rignot/rgadi'));
%% File directory and files content
file_dir = pwd;
km2m = 1000;
%work_dir = '/DFS-L/DATA/rignot/rgadi/08302022/run_13';
work_dir = pwd;
plot_dir = [work_dir filesep 'hFac_2D'];
if(~exist(plot_dir))
   mkdir(plot_dir);
end

delta_t = 1;
list_dyn = dir(fullfile(file_dir,'dynDiag*.data'));

day2s = 24 * 3600;
day2hr = 24;
nskip_hor = 2;
nskip_vert = 1;
scale_fac = 0.3;
fac_1 = 0.8;
x_max = 1000;
nfac = 3;
fac = 5;
skip = 1;
t_want = 10/24;
%% Co-ordinate values
Y = cumsum(squeeze(rdmds('DYC')))/km2m - 14.4;
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);


 for i=1:skip:length(list_dyn)
     code_T(i) = str2num(list_dyn(i).name(end-11:end-5)) * delta_t/day2s;
 end

 ind = find(min(abs(code_T - t_want))==abs(code_T - t_want));
%load('Max.mat');
figure('units','normalized','outerposition',[0 0 1 1]);
%% Making plots of 2D fields
for i=ind
%for i=length(list_T)
    %code_T = str2num(list_T(i).name(end-11:end-5)) * delta_t/day2s;
    code_T = t_want;
    data_2D = squeeze(rdmds(list_dyn(i).name(1:end-5)));
    hFac = squeeze(data_2D(:,:,3));
   
    hFac(hFac==0) = NaN;
    
    day = floor(code_T);
    hour = code_T - day;
    if(hour>0)
        time_vec = sprintf('%d Day %d hours',day,hour*day2hr);
    else
        time_vec = sprintf('%d Day',day);
    end
    
    clf;
    ax1 = gca;
    %ax1 = subplot(1,3,1);
    contourf(ax1,Ymesh',Zmesh',hFac,'LineStyle','none');
    %pcolorcen(Ymesh',Zmesh',T);
    colorbar;caxis(ax1,[0 1]);colormap(ax1,'jet');
    hold on;
    %quiver(ax1,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
    %   V(1:nskip_hor:end,1:nskip_vert:end)/(fac*Velmax),W(1:nskip_hor:end,1:nskip_vert:end)/(fac*Velmax),....
    %   scale_fac,'k','LineWidth',2,'MaxHeadSize',100);
    %axis equal;
    %hold off;
    %p = get(gca, 'Position');
    h = axes('Parent', gcf, 'OuterPosition', [0.1303 0.1605 0.1935 0.333]);
    contourf(h,Ymesh',Zmesh',hFac,'LineStyle','none');colorbar;caxis(h,[0 1]);colormap(h,'jet');
    %hold(h,'on');
    %quiver(h,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
    %    V(1:nskip_hor:end,1:nskip_vert:end),W(1:nskip_hor:end,1:nskip_vert:end),scale_fac*nfac,'k')
    %quiver(h,Ymesh',Zmesh',V,W,scale_fac*nfac,'k');
    %hold(h,'off');
    xlim(h,[-2 2]);
    ylim(h,[-480 -430]);
    %hold off;
    %xlim(h,[0 x_max]);
    %ylim(h,[-500 -400]);
    %p = get(h, 'Position');
    %h1 = axes('Parent', gcf, 'Position', [p(1)+0.05 p(2)+0.2 p(3)-0.25 p(4)-0.25]);
    %contourf(h1,Ymesh',Zmesh',T,'LineStyle','none');colorbar;caxis(h1,[Tmin Tmax]);colormap(h1,'jet');
    %hold(h1,'on');
    %quiver(h1,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
    %    V(1:nskip_hor:end,1:nskip_vert:end),W(1:nskip_hor:end,1:nskip_vert:end),scale_fac*nfac,'k')
    %quiver(h1,Ymesh',Zmesh',V,W,scale_fac*nfac,'k');
    %hold(h1,'off');
    %xlim(h1,[0 x_max*0.3]);
    %ylim(h1,[-500 -480]);
    xlabel(ax1,'Distance from GL (m)');
    ylabel(ax1,'Depth (m)');
    title(ax1,sprintf('Hfac at Time = %s',time_vec));
    set(ax1,'FontSize',22,'FontWeight','bold');
    print(gcf,'-dpng','-r300',fullfile(plot_dir,sprintf('Hfac_fields_%04d',i)));


    
end
