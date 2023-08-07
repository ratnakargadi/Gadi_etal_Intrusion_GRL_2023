%% This code is written for generating linear spacing co-ordinates. 
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));

%% Constants
km2m = 10^3;
deg2km = 111.4;
%% 
linear_spac_start = 18.88; % km from the GL where you want the linear spacing to start (in km)
linear_spac_end = 54.7; % km from the GL where you want the linear spacing to end (in km)
length_end = 63;
mesh_start = 0;
%% 
dx_start = 20; % the spacing at the starting point of linear spacing (in m)
dx_end = 500; % the spacing at the ending point of linear spacing (in m)

%% 
nx_start = floor(linear_spac_start * km2m/dx_start);
x = linspace(mesh_start,linear_spac_start,nx_start+1);

if(dx_start~=(x(end)-x(end-1))*km2m)
    dx_start = (x(end) - x(end-1))*km2m;
end

%% Linear Spacing algorithm
m = (dx_end - dx_start)/((linear_spac_end - linear_spac_start) * km2m -dx_end);
dx = dx_start;
%dx = 0;
xn(1) = 0;
count = 1;
while(xn(count)<=(linear_spac_end-linear_spac_start))
    xn(count+1) = dx/km2m + xn(count);
    %dx = dx + m * (xn(count+1)-xn(count)) * km2m;
    dx = dx * (m + 1);
    %dx = dx * m * km2m;
    gd(count) = dx;
    count = count + 1;
end
x = [x xn(2:end-1)+linear_spac_start];
    
%% Straight line segment
nx = floor((length_end - x(end))*km2m/dx_end);
x = [x(1:end-1) linspace(x(end),length_end,nx+1)];
delX = x/deg2km;
delXX = [diff(delX)];

namF='dely_spacing.bin';
%if kwr > 2,
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,hIce,prec);fclose(fid);
 writebin(namF,delXX);
 fprintf(' done\n');
%end
 disp(length(delXX))