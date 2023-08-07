clear all;
clc;
close all;

t_start = 0:22;
t_end = t_start + 1;

figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(t_start)
    ploter(t_start(i),t_end(i));
end