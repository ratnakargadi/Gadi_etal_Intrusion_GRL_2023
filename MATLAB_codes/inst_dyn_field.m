%% This code will plot Salinity, temperature and hFac values for the 
% instantaneous dynamic log files. This will help in seeing the way hfac
% changes.
clear all;
clc;
close all;

addpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab');
file_dir = pwd;
list = dir(fullfile(file_dir,'dynInst*.data'));

delta_t = 25;
ylims = [-225 -220];
salt_lims = [33 33+3];
temp_lims = [-3 3];

Y = squeeze(rdmds('YC'));
if(length(unique(Y))==1)
    Y = squeeze(rdmds('YC'));
end

Z = -cumsum(squeeze(rdmds('DRF')));
[Ymesh,Zmesh] = meshgrid(Y,Z);

figure('units','normalized','outerposition',[0 0 1 1]);

for i=1:length(list)
    Data = squeeze(rdmds(list(i).name(1:end-5)));
    salt = squeeze(Data(:,:,5));
    temp = squeeze(Data(:,:,4));
    hFacS = squeeze(Data(:,:,end));
    salt(hFacS==0) = NaN;
    temp(hFacS==0) = NaN;
    disp(sprintf('Timestep %d extrem salt: [%4.2f %4.2f]',i,min(salt(:)),max(salt(:))));
    salt_extr(i,:) = [min(salt(:)) max(salt(:))];
    if(~isempty(salt_lims))
        salt(find(salt<salt_lims(1)|salt>salt_lims(2))) = NaN;
    end
    
    if(~isempty(temp_lims))
        temp(find(temp<temp_lims(1)|temp>temp_lims(2))) = NaN;
    end
    
    clf;
    ax1 = subplot(1,3,1);
    contourf(ax1,Ymesh,Zmesh,salt');
    xlabel(ax1,'Longitude (\circ)');
    ylabel(ax1,'Depth (m)');
    ylim(ax1,[ylims]);
    colorbar;
    caxis(ax1,[salt_lims]);
    title(ax1,'Salinity');
    
    ax2 = subplot(1,3,2);
    contourf(ax2,Ymesh,Zmesh,temp');
    xlabel(ax2,'Longitude (\circ)');
    ylabel(ax2,'Depth (m)');
    ylim(ax2,[ylims]);
    colorbar;
    caxis(ax2,[temp_lims]);
    title(ax2,'Temperature');
    
    ax3 = subplot(1,3,3);
    contourf(ax3,Ymesh,Zmesh,hFacS');
    xlabel(ax3,'Longitude (\circ)');
    ylabel(ax3,'Depth (m)');
    ylim(ax3,[ylims]);
    colorbar;
    caxis(ax3,[0 1]);
    title(ax3,'hFac');
    
    suptitle(sprintf('Instantaneous dynamic field at t=%d s',delta_t*str2num(list(i).name(end-8:end-5))));
    str = ['dynInst' '_' 'timestamp' '_' list(i).name(end-8:end-5)];
    print(gcf,'-dpng','-r300',fullfile(pwd,str));
end
