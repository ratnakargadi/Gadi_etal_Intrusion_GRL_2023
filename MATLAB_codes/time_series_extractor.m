%% This code extracts the time-series of salinity and temperature at maximum
% 8 locations specified. This code is used to check for the convergence of
% the oceanic fields in the ice-shelf cavities
% COMMENT THE NEXT THREE LINES IF USING A SHELL SCRIPT
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%% FILE DIRECTORY AND LOCATION OF FILES
file_dir = pwd;
delta_t = 4;
day2s = 24 * 3600;
ver = 2;
list_T = dir(fullfile(file_dir,'T*.data'));
list_S = dir(fullfile(file_dir,'S.0*data'));
list_V = dir(fullfile(file_dir,'V*.data'));
list_U = dir(fullfile(file_dir,'U*.data'));
list_W = dir(fullfile(file_dir,'W*.data'));
n_time = length(list_T);
topofile = [file_dir filesep 'shelficeTopo.Lin.bin'];
for i=1:length(list_S)
    code_t(i) = str2num(list_S(i).name(end-10:end-5));
end
code_t = code_t * delta_t/day2s;

%% Reading the co-ordinates of the ice-shelf file
Y = cumsum(squeeze(rdmds('DYC')));
Z = -cumsum(squeeze(rdmds('DRF')));
dZ = max(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);
Ymesh = Ymesh';Zmesh = Zmesh';

%% Y-Z (X-Z) locations where time-series needs to be extracted
X_loc = [100 500 1000 1500 1500 1000 500 100];
H = readbin(topofile,[1 length(Y)])';
Z_loc = [ones(1,4)*min(H) interp1(Y,H,X_loc(5:end))-1.5 * dZ];
Z_loc(end) = min(H);
time_vec.t = zeros(n_time,length(X_loc));
time_vec.s = zeros(size(time_vec.t));
time_vec.u = zeros(size(time_vec.t));
time_vec.v = zeros(size(time_vec.t));
time_vec.w = zeros(size(time_vec.t));

func_interp = @(x,y,field)scatteredInterpolant(x(:),y(:),field(:));
%% Running a loop to extract the time-series
for i=1:n_time
    T = squeeze(rdmds(list_T(i).name(1:end-5)));
    S = squeeze(rdmds(list_S(i).name(1:end-5)));
    U = squeeze(rdmds(list_U(i).name(1:end-5)));
    V = squeeze(rdmds(list_V(i).name(1:end-5)));
    W = squeeze(rdmds(list_W(i).name(1:end-5)));
    
    S(S==0) = NaN;
    T(isnan(S)) = NaN;
    U(isnan(S)) = NaN;
    V(isnan(S)) = NaN;
    W(isnan(S)) = NaN;
    
    func = func_interp(Ymesh,Zmesh,T);
    time_vec.t(i,:) = func(X_loc,Z_loc);
    clear func;
    
    func = func_interp(Ymesh,Zmesh,S);
    time_vec.s(i,:) = func(X_loc,Z_loc);
    clear func;
    
    func = func_interp(Ymesh,Zmesh,U);
    time_vec.u(i,:) = func(X_loc,Z_loc);
    clear func;
    
    func = func_interp(Ymesh,Zmesh,V);
    time_vec.v(i,:) = func(X_loc,Z_loc);
    clear func;
    
    func = func_interp(Ymesh,Zmesh,W);
    time_vec.w(i,:) = func(X_loc,Z_loc);
    clear func;
    
    disp(sprintf('%d out of %d time steps',i,n_time));
end

%% Plotting the time-series extracted
figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(X_loc)
    ax = subplot(4,2,i);
    plot(ax,code_t,time_vec.t(:,i));
    xlabel('Time (days)');
    ylabel('Temp (C)');
    title(sprintf('At (Y,Z) = (%d m,%4.2f m)',X_loc(i),Z_loc(i)));
end
print(gcf,'-dpng','-r300',fullfile(file_dir,sprintf('Temp_time_series_%d',ver)));

figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(X_loc)
    ax = subplot(4,2,i);
    plot(ax,code_t,time_vec.s(:,i));
    xlabel('Time (days)');
    ylabel('Salt (psu)');
    title(sprintf('At (Y,Z) = (%d m,%4.2f m)',X_loc(i),Z_loc(i)));
end
print(gcf,'-dpng','-r300',fullfile(file_dir,sprintf('Salt_time_series_%d',ver)));

figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(X_loc)
    ax = subplot(4,2,i);
    plot(ax,code_t,time_vec.u(:,i));
    xlabel('Time (days)');
    ylabel('U (m/s)');
    title(sprintf('At (Y,Z) = (%d m,%4.2f m)',X_loc(i),Z_loc(i)));
end
print(gcf,'-dpng','-r300',fullfile(file_dir,sprintf('Zonal_time_series_%d',ver)));

figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(X_loc)
    ax = subplot(4,2,i);
    plot(ax,code_t,time_vec.v(:,i));
    xlabel('Time (days)');
    ylabel('V (m/s)');
    title(sprintf('At (Y,Z) = (%d m,%4.2f m)',X_loc(i),Z_loc(i)));
end
print(gcf,'-dpng','-r300',fullfile(file_dir,sprintf('Meridonal_time_series_%d',ver)));

figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(X_loc)
    ax = subplot(4,2,i);
    plot(ax,code_t,time_vec.w(:,i));
    xlabel('Time (days)');
    ylabel('W (m/s)');
    title(sprintf('At (Y,Z) = (%d m,%4.2f m)',X_loc(i),Z_loc(i)));
end
print(gcf,'-dpng','-r300',fullfile(file_dir,sprintf('Vertical_time_series_%d',ver)));
