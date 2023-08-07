function [Tb,Sb,wb] = meltrates(Tm,Sm,topofile,gammaS,gammaT)
% This function computes the melt rate under the ice sheet using the
% coefficients as specified in "Modeling Thermodynamic Ice-Ocean Interactions 
% at the Base of an Ice Shelf" (Holland and Jenkins 1999)
nx = size(squeeze(rdmds('DXC')));
rhoref = 1025;
rhoI = 920;
a = -5.73 * 10^-2;
b = 9.39 * 10^-2;
c = -7.53 * 10^-8;
cpI = 2009;
cpM = 3974;
g = 9.81;
Lf = 3.34 * 10^5;
Kit = 1.14 * 10^-6;
Kmt = 1;
Ts = -25;
Si = 0;
y2s = 365 * 24 * 3600;
prec = 'real*8';
fileId = fopen(topofile,'r','b');
H = abs(fread(fileId,nx,prec));
fclose(fileId);
Pb = H * rhoref * g;
hFac = squeeze(rdmds('hFacC'));
for i=1:size(hFac,1)
    [row,col] = find(hFac(i,:)>0,1,'first');
    if(~isempty(col))
        hval(i) = hFac(i,col);
    else
        hval(i) = 1;
    end
end
dz = squeeze(rdmds('DRC'));
dz = max(dz);
h = hval * dz';
adash=(rhoI*cpI*Kit./H-rhoref*gammaT*cpM)*a;
bdash=(b+c*Pb-Ts)./H*rhoI*cpI*Kit;
bdash=-(b+c*Pb-Tm)*gammaT*rhoref*cpM+bdash;
bdash=bdash+rhoref*gammaS*Lf;
cdash=rhoref*gammaS*Sm*Lf;
Sb=(-bdash+sqrt(bdash.*bdash+4*adash.*cdash))./(2*adash);
Tb = a * Sb + b + c * Pb;
wb = rhoref * gammaS * (Sb - Sm)./(rhoI * (Si - Sb)) * y2s;
end


