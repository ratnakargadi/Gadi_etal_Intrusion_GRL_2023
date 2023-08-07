%% This code plots the instantaneous velocity along the GL
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
%%
se = squeeze(rdmds('surfDiag',NaN));
GL = squeeze(se(:,7,:));
eta = squeeze(se(:,1,:));
GL_inst = GL(1,:);
file_GL = dir(fullfile(pwd,'surfDiag*.data'));
file_s = dir(fullfile(pwd,'V*.data'));
deltaT = 1;
sec2hr = 3600;

for i=1:length(file_GL)
    t_GL(i) = str2num(file_GL(i).name(10:19)) * deltaT/sec2hr;
end

for i=1:length(file_s)
    t_V(i) = str2num(file_s(i).name(3:12)) * deltaT/sec2hr;
end

GL_v = interp1(t_GL,GL_inst,t_V);
%% 
V = squeeze(rdmds('V',NaN));
W = squeeze(rdmds('W',NaN));
YC = squeeze(rdmds('YC'));
DYC = squeeze(rdmds('DYC'));
y_ind = find(YC<=GL_inst(1));
Z = 1:size(V,2);

[Ymesh,Zmesh] = meshgrid(YC,Z);
[Ywant,Zwant] = meshgrid(YC(y_ind),Z);
for i=1:size(V,3)
    
    Dum = interp2(Ymesh,Zmesh,squeeze(V(:,:,i))',YC(i),Z);
    for j=1:size(Dum,1)
        D(j) = sum(Dum(j,:) .* DYC(y_ind))./sum(DYC(y_ind));
    end
    DW = interp2(Ymesh,Zmesh,squeeze(W(:,:,i))',GL_v(i),Z);
    try
        Vvel(i) = D(min(find(abs(D)>0)));
    catch
        Vvel(i) = 0;
    end
    try
        Wvel(i) = DW(min(find(abs(DW)>0)));
    catch
        Wvel(i) = 0;
    end
end

t_indd = find(t_V>=0&t_V<=6);
int_Vvel = sum(Vvel(t_indd))/length(t_indd);
disp(int_Vvel);
t_ind = find(t_V>=6&t_V<=12);
int_Vvel =sum(Vvel(t_ind))/length(t_ind);
disp(int_Vvel);
figure('units','normalized','outerposition',[0 0 1 1]);
plot(t_V,Vvel,'k','LineWidth',2);
hold on;
plot(t_V,Wvel,'g','LineWidth',2);
hold off;
xlabel('Time (hr)');
ylabel('Velocity (m/s)');
legend('Horizontal Vel.','Vertical Vel');
title('Velocities at instantaneous GL');
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'inst_vel'));

figure;
plot(t_GL,eta(1190,:));
xlabel('Time (hr)');
ylabel('\eta at domain end (m)');
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'inst_eta'));

