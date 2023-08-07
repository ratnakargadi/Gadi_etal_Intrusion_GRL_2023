%% This code will plot the 2D velocity fields time series under an ice-sheet
% Comment the next three lines after debugging
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
addpath(genpath('/DFS-L/DATA/rignot/rgadi'));
%% File directory and files content
file_dir = pwd;
plot_dir = [file_dir filesep 'twoD_fields'];
if(~exist(plot_dir))
   mkdir(plot_dir);
end

delta_t = 3;
list_T = dir(fullfile(file_dir,'T*.data'));
list_S = dir(fullfile(file_dir,'S.0*.data'));
list_U = dir(fullfile(file_dir,'U*.data'));
list_V = dir(fullfile(file_dir,'V*.data'));
list_W = dir(fullfile(file_dir,'W*.data'));
day2s = 24 * 3600;
day2hr = 24;
nskip_hor = 3;
nskip_vert = 10;
scale_fac = 0.3;
x_max = 1000;
nfac = 3;
%% Co-ordinate values
Y = cumsum(squeeze(rdmds('DYC')));
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);

%Umax = -1000;Umin = 1000;
Vmax = -1000;Vmin = 1000;
Wmax = -1000;Wmin = 1000;
Velmin = 1000;Velmax = -1000;
Smin = 1000;Smax = -1000;
Tmin = 1000;Tmax = -1000;

for i=1:length(list_T)
    S = squeeze(rdmds(list_S(i).name(1:end-5)));
    S(S==0) = NaN;
    S(end,:) = NaN;
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
    if(mod(i,100)==0)
        disp(i);
    end
end

Vel_min = Velmin;
Vel_max = Velmax;

save('Max.mat','Vel_min','Vel_max','Smin','Smax','Tmin','Tmax','Wmin','Wmax','Vmin','Vmax');