%% 
clear all;
clc;
close all;

%%
km2m = 10^3;
file_1 = [pwd filesep 'tidal_average_5_1.mat'];
file_2 = [pwd filesep 'tidal_average_5_2.mat'];
file_3 = [pwd filesep 'tidal_average_5_3.mat'];
file_4 = [pwd filesep 'tidal_average_5_4.mat'];
file_5 = [pwd filesep 'tidal_average_5_5.mat'];
file_6 = [pwd filesep 'tidal_average_5_6.mat'];

%%
A_1 = load(file_1);
A_2 = load(file_2);
A_3 = load(file_3);
A_4 = load(file_4);
A_5 = load(file_5);
A_6 = load(file_6);

ind = find(A_1.Y<0);
A_1.int_melt = sum(A_1.Y(ind).* A_1.melt_tidal_avg(ind))/sum(A_1.Y(ind)) * 1.028;
A_2.int_melt = sum(A_2.Y(ind).* A_2.melt_tidal_avg(ind))/sum(A_2.Y(ind)) * 1.028;
A_3.int_melt = sum(A_3.Y(ind).* A_3.melt_tidal_avg(ind))/sum(A_3.Y(ind)) * 1.028;
A_4.int_melt = sum(A_4.Y(ind).* A_4.melt_tidal_avg(ind))/sum(A_4.Y(ind)) * 1.028;
A_5.int_melt = sum(A_5.Y(ind).* A_5.melt_tidal_avg(ind))/sum(A_5.Y(ind)) * 1.028;
A_6.int_melt = sum(A_6.Y(ind).* A_6.melt_tidal_avg(ind))/sum(A_6.Y(ind)) * 1.028;


%%
figure('units','normalized','OuterPosition',[0 0 1 1]);
plot(A_1.Y/km2m,A_1.melt_tidal_avg,'r','LineWidth',2);
hold on;
plot(A_2.Y/km2m,A_2.melt_tidal_avg,'k','LineWidth',2);
hold on;
plot(A_3.Y/km2m,A_3.melt_tidal_avg,'b','LineWidth',2);
hold on;
plot(A_4.Y/km2m,A_4.melt_tidal_avg,'m','LineWidth',2);
hold on;
plot(A_5.Y/km2m,A_5.melt_tidal_avg,'c','LineWidth',2);
hold on;
plot(A_6.Y/km2m,A_6.melt_tidal_avg,'Color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on;
plot([0 0],[0 60],'y','LineWidth',3,'LineStyle','--');
hold on;
plot([-5 -5],[0 60],'g','LineWidth',3,'LineStyle','--');
hold off;
L = legend('1st cycle','2nd cycle','3rd cycle','4th cycle','5th cycle','6th cycle','Low tide','High tide');
L.FontSize = 22;
L.FontWeight = 'bold';
xlabel('Distance from GL (km)');
ylabel('Basal Melt rate (m/yr)');
title('Tidally averaged basal melt rate');
set(gca,'FontSize',22,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'tidal_average_cycle_plot'));

figure;
plot([1:6],[A_1.int_melt A_2.int_melt A_3.int_melt A_4.int_melt A_5.int_melt A_6.int_melt],'k','LineWidth',2);
hold on;
plot([1:6],[max(A_1.melt_tidal_avg) max(A_2.melt_tidal_avg) max(A_3.melt_tidal_avg) max(A_4.melt_tidal_avg) max(A_5.melt_tidal_avg) max(A_6.melt_tidal_avg)],'r','LineWidth',2);
hold off;
xlabel('Tidal cycle');
ylabel('Melt rate (m/yr)');
L = legend('GZ avg.','Max');
L.FontSize = 16;
L.FontWeight = 'bold';
set(gca,'FontSize',16,'FontWeight','bold');
print(gcf,'-dpng','-r300',fullfile(pwd,'GZ_convergence'));
