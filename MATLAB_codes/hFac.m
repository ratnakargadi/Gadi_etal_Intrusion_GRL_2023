clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));

data = squeeze(rdmds('dynDiag',NaN));
hFac = squeeze(data(:,:,3,:));
hFac(hFac==0) = NaN;
km2m = 10^3;
GZ = 6;
Y = cumsum(squeeze(rdmds('DYC')))/km2m - GZ;
loc_want = [-GZ+1 -GZ/2 0];
delta_t = 1;
sec2hr = 3600;
list = dir(fullfile(pwd,'dynDiag.*data'));

for k=1:length(list)
    code_T(k) = str2num(list(k).name(end-11:end-5)) * delta_t/sec2hr;
end

for j=1:length(loc_want)
    ind(j) = min(find(min(abs(Y - loc_want(j)))==abs(Y - loc_want(j))));
end

for i=1:size(hFac,3)
    hFac_2d = squeeze(hFac(:,:,i));
    for j=1:length(ind)
        ind_nan(i,j) = hFac_2d(ind(j),min(find(~isnan(hFac_2d(ind(j),:)))));
    end
end


figure('units','normalized','outerposition',[0 0 1 1]);
plot(code_T,ind_nan(:,1),'b','LineWidth',2)
hold on
plot(code_T,ind_nan(:,2),'r','LineWidth',2)
hold on
plot(code_T,ind_nan(:,3),'k','LineWidth',2)
L = legend('GL - 5 km','GL - 3 km','GL');
L.FontSize = 16;
L.FontWeight = 'bold';
ylim([0 1]*max(abs(ind_nan(:)))*1.2);
L.Location = 'best';
xlabel('Time (hrs)');
ylabel('Mom (\rho*V)');
xlim([0 100]);
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'hFac_1'));
