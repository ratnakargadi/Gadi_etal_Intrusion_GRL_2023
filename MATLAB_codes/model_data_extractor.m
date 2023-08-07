%% This code will find the recent file and extract the datasets at the requested 
% locations
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%% File location and directory
file_dir = pwd;
file_save_dir = [file_dir filesep 'Model_data'];
if(~exist(file_save_dir))
    mkdir(file_save_dir);
end

list_S = dir(fullfile(file_dir,'S.0*'));
list_T = dir(fullfile(file_dir,'T.0*'));
S = squeeze(rdmds('S',str2double(list_S(end).name(end-12:end-5))));
S(S==0) = NaN;
T = squeeze(rdmds('T',str2double(list_T(end).name(end-12:end-5))));
T(isnan(S)) = NaN;
km2m = 1000;

%% Co-ordinate values
Y = cumsum(squeeze(rdmds('DYC')))/km2m;
Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);

%% Locations where you want to extract
y_want = [3;16;26];
z_want = Z;
[Y_W,Z_W] = meshgrid(y_want,z_want);

if(size(S,1)~=size(Ymesh,1))
    S  = S'; T = T';
end

func_S = scatteredInterpolant(Ymesh(:),Zmesh(:),S(:));
S_w = func_S(Y_W,Z_W);

func_T = scatteredInterpolant(Ymesh(:),Zmesh(:),T(:));
T_w = func_T(Y_W,Z_W);

Pressure = abs(Z);

for i=1:length(y_want)
    filename_save = [file_save_dir filesep sprintf('Model_%02s.mat',num2str(y_want(i),'%02d'))];
    Salinity = S_w(:,i);
    Temperature = T_w(:,i);
    
    save(filename_save,'Salinity','Temperature','Pressure');
    clear Salinity Temperature;
end