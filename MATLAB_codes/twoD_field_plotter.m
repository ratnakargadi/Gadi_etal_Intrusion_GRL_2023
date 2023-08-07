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
plot_dir = [work_dir filesep 'twoD_fields'];
if(~exist(plot_dir))
   mkdir(plot_dir);
end

delta_t = 1;
list_T = dir(fullfile(file_dir,'T*.data'));
list_S = dir(fullfile(file_dir,'S.0*.data'));
list_U = dir(fullfile(file_dir,'U*.data'));
list_V = dir(fullfile(file_dir,'V*.data'));
list_W = dir(fullfile(file_dir,'W*.data'));
day2s = 24 * 3600;
day2hr = 24;
nskip_hor = 2;
nskip_vert = 1;
scale_fac = 0.3;
fac_1 = 0.8;
x_max = 1000;
nfac = 3;
fac = 5;
skip = 5;
%% Co-ordinate values
Y = cumsum(squeeze(rdmds('DYC')))/km2m;
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);

%Umax = -1000;Umin = 1000;
Vmax = -1000;Vmin = 1000;
Wmax = -1000;Wmin = 1000;
Velmin = 1000;Velmax = -1000;
Smin = 1000;Smax = -1000;
Tmin = 1000;Tmax = -1000;

 for i=1:1:length(list_T)
     S = squeeze(rdmds(list_S(i).name(1:end-5)));
     S(S==0) = NaN;
     %S(end,:) = NaN;
     T = squeeze(rdmds(list_T(i).name(1:end-5)));
     T(isnan(S)) = NaN;
     %U = squeeze(rdmds(list_U(i).name(1:end-5)));
     V = squeeze(rdmds(list_V(i).name(1:end-5)));
     W = squeeze(rdmds(list_W(i).name(1:end-5)));
     Vel = sqrt(V.^2 + W.^2);
     
     %U(isnan(S)) = NaN;
     V(isnan(S)) = NaN;
     W(isnan(S)) = NaN;
     
     %Umax = max(Umax,max(U(:)));
     %Umin = min(Umin,min(U(:)));
     Vmax = max(Vmax,max(V(:)));
     Vmin = min(Vmin,min(V(:)));
     Wmax = max(Wmax,max(W(:)));
     Wmin = min(Wmin,min(W(:)));
     Velmin = min(Velmin,min(Vel(:)));
     Velmax = max(Velmax,max(Vel(:)));
     Tmin = min(Tmin,min(T(:)));
     Tmax = max(Tmax,max(T(:)));
     Smin = min(Smin,min(S(:)));
     Smax = max(Smax,max(S(:)));
 end
 
 Vel_min = Velmin;
 Vel_max = Velmax;
%load('Max.mat');
figure('units','normalized','outerposition',[0 0 1 1]);
%% Making plots of 2D fields
for i=1:skip:length(list_T)
%for i=length(list_T)
    code_T = str2num(list_T(i).name(end-11:end-5)) * delta_t/day2s;
    T = squeeze(rdmds(list_T(i).name(1:end-5)));
    S = squeeze(rdmds(list_S(i).name(1:end-5)));
    %U = squeeze(rdmds(list_U(i).name(1:end-5)));
    V = squeeze(rdmds(list_V(i).name(1:end-5)));
    W = squeeze(rdmds(list_W(i).name(1:end-5)));
    
    S(S==0) = NaN;
    S(end,:) = NaN;
    T(isnan(S)) = NaN;
    V(isnan(S)) = NaN;
    %U(isnan(S)) = NaN;
    W(isnan(S)) = NaN;
    
    day = floor(code_T);
    hour = code_T - day;
    if(hour>0)
        time_vec = sprintf('%d Day %d hours',day,hour*day2hr);
    else
        time_vec = sprintf('%d Day',day);
    end
    
    clf;
%     ax1 = gca;
%     %ax1 = subplot(1,3,1);
%     contourf(ax1,Ymesh',Zmesh',T,'LineStyle','none');
%     %pcolorcen(Ymesh',Zmesh',T);
%     colorbar;caxis(ax1,[Tmin Tmax]);colormap(ax1,'jet');
%     hold on;
%     quiver(ax1,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
%        V(1:nskip_hor:end,1:nskip_vert:end)/(fac*Velmax),W(1:nskip_hor:end,1:nskip_vert:end)/(fac*Velmax),....
%        scale_fac,'k','LineWidth',2,'MaxHeadSize',100);
%     %axis equal;
%     hold off;
%     %p = get(gca, 'Position');
%     %h = axes('Parent', gcf, 'Position', [p(1)+0.05 p(2)+0.4 p(3)-0.25 p(4)-0.45]);
%     %contourf(h,Ymesh',Zmesh',T,'LineStyle','none');colorbar;caxis(h,[Tmin Tmax]);colormap(h,'jet');
%     %hold(h,'on');
%     %quiver(h,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
%     %    V(1:nskip_hor:end,1:nskip_vert:end),W(1:nskip_hor:end,1:nskip_vert:end),scale_fac*nfac,'k')
%     %quiver(h,Ymesh',Zmesh',V,W,scale_fac*nfac,'k');
%     %hold(h,'off');
%     %xlim(h,[0 x_max]);
%     %ylim(h,[-500 -400]);
%     %p = get(h, 'Position');
%     %h1 = axes('Parent', gcf, 'Position', [p(1)+0.05 p(2)+0.2 p(3)-0.25 p(4)-0.25]);
%     %contourf(h1,Ymesh',Zmesh',T,'LineStyle','none');colorbar;caxis(h1,[Tmin Tmax]);colormap(h1,'jet');
%     %hold(h1,'on');
%     %quiver(h1,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
%     %    V(1:nskip_hor:end,1:nskip_vert:end),W(1:nskip_hor:end,1:nskip_vert:end),scale_fac*nfac,'k')
%     %quiver(h1,Ymesh',Zmesh',V,W,scale_fac*nfac,'k');
%     %hold(h1,'off');
%     %xlim(h1,[0 x_max*0.3]);
%     %ylim(h1,[-500 -480]);
%     xlabel(ax1,'Distance from GL (m)');
%     ylabel(ax1,'Depth (m)');
%     title(ax1,sprintf('Temp (C) at Time = %s',time_vec));
%     print(gcf,'-dpng','-r300',fullfile(plot_dir,sprintf('Temp_fields_%04d',i)));
%     
%     clf;
%     ax2 = gca;
%     %ax2 = subplot(1,3,2);
%     contourf(ax2,Ymesh',Zmesh',S,'LineStyle','none');
%     %pcolorcen(Ymesh',Zmesh',S);
%     colorbar;caxis(ax2,[Smin Smax]);colormap(ax2,'jet');
%     hold on;
%     quiver(ax2,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
%        V(1:nskip_hor:end,1:nskip_vert:end)/(fac*Velmax),W(1:nskip_hor:end,1:nskip_vert:end)/(fac*Velmax),....
%        scale_fac,'k','LineWidth',2,'MaxHeadSize',100);
%     hold off;
%     %p = get(gca, 'Position');
%     %h = axes('Parent', gcf, 'Position', [p(1)+0.05 p(2)+0.4 p(3)-0.25 p(4)-0.45]);
%     %contourf(h,Ymesh',Zmesh',S,'LineStyle','none');colorbar;caxis(h,[Smin Smax]);colormap(h,'jet');
%     %hold(h,'on');
%     %quiver(h,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
%     %    V(1:nskip_hor:end,1:nskip_vert:end),W(1:nskip_hor:end,1:nskip_vert:end),scale_fac*nfac,'k')
%     %quiver(h,Ymesh',Zmesh',V,W,scale_fac*nfac,'k');
%     %hold(h,'off');
%     %xlim(h,[0 x_max]);
%     %ylim(h,[-500 -400]);
%     %p = get(h, 'Position');
%     %h1 = axes('Parent', gcf, 'Position', [p(1)+0.05 p(2)+0.2 p(3)-0.25 p(4)-0.25]);
%     %contourf(h1,Ymesh',Zmesh',S,'LineStyle','none');colorbar;caxis(h1,[Smin Smax]);colormap(h1,'jet');
%     %hold(h1,'on');
%     %quiver(h1,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
%     %    V(1:nskip_hor:end,1:nskip_vert:end),W(1:nskip_hor:end,1:nskip_vert:end),scale_fac*nfac,'k')
%     %quiver(h1,Ymesh',Zmesh',V,W,scale_fac*nfac,'k');
%     %hold(h1,'off');
%     %xlim(h1,[0 x_max*0.3]);
%     %ylim(h1,[-500 -480]);
%     xlabel(ax2,'Distance from GL (m)');
%     ylabel(ax2,'Depth (m)');
%     title(ax2,sprintf('Salt (psu) at Time = %s',time_vec));
%     print(gcf,'-dpng','-r300',fullfile(plot_dir,sprintf('Salt_fields_%04d',i)));
    
    clf;
    ax3 = gca;
    %ax3 = subplot(1,3,3);
    %contourf(ax3,Ymesh',Zmesh',V,'LineStyle','none');
    pcolorcen(Ymesh',Zmesh',V);
    colorbar;caxis(ax3,[-1 1]*max(abs(Vmin),abs(Vmax))*fac_1);colormap(ax3,cmocean('balance','pivot',0));
    hold on;
    quiver(ax3,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
       V(1:nskip_hor:end,1:nskip_vert:end)/(fac*Velmax),W(1:nskip_hor:end,1:nskip_vert:end)/(fac*Velmax),....
       scale_fac,'k','LineWidth',2,'MaxHeadSize',100);
    hold off;
    %p = get(gca, 'Position');
    %h = axes('Parent', gcf, 'Position', [p(1)+0.05 p(2)+0.4 p(3)-0.25 p(4)-0.45]);
    %contourf(h,Ymesh',Zmesh',V,'LineStyle','none');colorbar;caxis(h,[-1 1]*max(abs(Vmin),abs(Vmax)));colormap(h,cmocean('balance','pivot',0));
    %hold(h,'on');
    %quiver(h,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
    %    V(1:nskip_hor:end,1:nskip_vert:end),W(1:nskip_hor:end,1:nskip_vert:end),scale_fac*nfac,'k')
    %quiver(h,Ymesh',Zmesh',V,W,scale_fac*nfac,'k');
    %hold(h,'off');
    %xlim(h,[0 x_max]);
    %ylim(h,[-500 -400]);
    %p = get(h, 'Position');
    %h1 = axes('Parent', gcf, 'Position', [p(1)+0.05 p(2)+0.2 p(3)-0.25 p(4)-0.25]);
    %contourf(h1,Ymesh',Zmesh',V,'LineStyle','none');colorbar;caxis(h1,[-1 1]*max(abs(Vmin),abs(Vmax)));colormap(h1,cmocean('balance','pivot',0));
    %hold(h1,'on');
    %quiver(h1,Ymesh(1:nskip_vert:end,1:nskip_hor:end)',Zmesh(1:nskip_vert:end,1:nskip_hor:end)',....
    %    V(1:nskip_hor:end,1:nskip_vert:end),W(1:nskip_hor:end,1:nskip_vert:end),scale_fac*nfac,'k')
    %quiver(h1,Ymesh',Zmesh',V,W,scale_fac*nfac,'k');
    %hold(h1,'off');
    %xlim(h1,[0 x_max*0.3]);
    %ylim(h1,[-500 -480]);
    xlabel(ax3,'Distance from GL (m)');
    ylabel(ax3,'Depth (m)');
    title(ax3,sprintf('2D Vel (m/s) at Time = %s',time_vec));
    print(gcf,'-dpng','-r300',fullfile(plot_dir,sprintf('Vel_fields_%04d',i)));
    
%     ax3 = subplot(2,3,3);
%     contourf(ax3,Ymesh',Zmesh',U,'LineStyle','none');colorbar;caxis(ax3,[Umin Umax]);
%     xlabel(ax3,'Distance from GL (m)');
%     ylabel(ax3,'Depth (m)');
%     title(ax3,'Zonal velocity (m/s)');
%     
%     ax4 = subplot(2,3,4);
%     contourf(ax4,Ymesh',Zmesh',V,'LineStyle','none');colorbar;caxis(ax4,[Vmin Vmax]);
%     xlabel(ax4,'Distance from GL (m)');
%     ylabel(ax4,'Depth (m)');
%     title(ax4,'Meridonal velocity (m/s)');
%     
%     ax5 = subplot(2,3,[5]);
%     contourf(ax5,Ymesh',Zmesh',W,'LineStyle','none');colorbar;caxis(ax5,[Wmin Wmax]);
%     xlabel(ax5,'Distance from GL (m)');
%     ylabel(ax5,'Depth (m)');
%     title(ax5,'Vertical velocity');
    
    
end
