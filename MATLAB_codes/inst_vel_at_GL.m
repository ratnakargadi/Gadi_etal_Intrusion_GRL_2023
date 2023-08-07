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
Z = 1:size(V,2);

[Ymesh,Zmesh] = meshgrid(YC,Z);
for i=1:size(V,3)    
    D = interp2(Ymesh,Zmesh,squeeze(V(:,:,i))',GL_v(i),Z);
    DD = interp2(Ymesh,Zmesh,squeeze(V(:,:,i))',GL_inst(1),Z);
    DDD = interp2(Ymesh,Zmesh,squeeze(V(:,:,i))',(max(GL_v)+min(GL_v))/2,Z);
    DDDD = interp2(Ymesh,Zmesh,squeeze(V(:,:,i))',min(GL_v),Z);
    DW = interp2(Ymesh,Zmesh,squeeze(W(:,:,i))',GL_v(i),Z);
    try
        Vvel(i) = D(min(find(abs(D)>0)));
    catch
        Vvel(i) = 0;
    end
    try
        VVvel(i) = DD(min(find(abs(DD)>0)));
    catch
        VVvel(i) = 0;
    end
    try
        VVVvel(i) = DDD(min(find(abs(DDD)>0)));
    catch
        VVVvel(i) = 0;
    end
    try
        VVVVvel(i) = DDDD(min(find(abs(DDDD)>0)));
    catch
        VVVVvel(i) = 0;
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
plot(t_V,VVvel,'g','LineWidth',2);
hold on;
plot(t_V,VVVvel,'r','LineWidth',2);
hold on;
plot(t_V,VVVVvel,'b','LineWidth',2);
hold off;
xlabel('Time (hr)');
ylabel('Horizontal Velocity (m/s)');
legend('Inst. GL','Cavity opening','half cavity length','Cavity length');
title('Velocities at instantaneous GL');
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'vel_2D'));

figure;
plot(t_GL,eta(1190,:));
xlabel('Time (hr)');
ylabel('\eta at domain end (m)');
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'inst_eta'));

