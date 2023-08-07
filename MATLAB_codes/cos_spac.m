function F = cos_spac(x,lx_1,lx_2)
%% This function is optimized basing on the choices of beginning mesh size and
% ending mesh size
F1 = [x(1)/2 * (1 - cos(pi/(2 * x(2)))) - lx_1 ....
    x(1)/2 * (1 - cos((x(2) - 1) * pi/(2 * x(2)))) - lx_2];
F = sqrt(sum(F1.^2)/length(F1));
end

