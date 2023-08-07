% This is a matlab script that generates the input data
% require matlab functions for equation of state
clear all;
clc;
close all;

addpath(genpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab'));
kwr=3;
% Dimensions of grid
nx=1;
ny=1200;
nr=2400;
delz=10;
km2m = 10^3;
deg2km = 111.4;
calving_front_ice_side = -6;
opt_topo = 1; % Use 0 for simplified triangular geometry
opt_ini = 1; % Use 0 for simple uniform initialization
%data_smooth = [pwd filesep 'Smoothed_data_7_shifted_36.mat']; % provide it when using real time data
data_smooth = [pwd filesep 'petermann_modified_6.mat'];
if(~isempty(data_smooth))
    load(data_smooth);
    clear bed_smooth;
    base = load('petermann_modified.mat','dist2GL','bed_smooth');
    bed_smooth = interp1(base.dist2GL,base.bed_smooth,dist2GL);
%     if(exist('GL_pos','var'))
%         %dist2GL = dist2GL - GL_pos;
%         bed_GL = interp1(dist2GL,bed_smooth,0);
%         base_GL = interp1(dist2GL,base_smooth,0);
%         bed_smooth(dist2GL<0) = [];
%         base_smooth(dist2GL<0) = [];
%         dist2GL(dist2GL<0) = [];
%         bed_smooth = [bed_GL;bed_smooth];
%         base_smooth = [base_GL;base_smooth];
%         dist2GL = [0;dist2GL];
%     end
%     load('mesh_spacing.mat');
%     if(length(meshh)<ny)
%         meshh = [meshh ones(1,ny-length(meshh))*meshh(end)];
%     elseif(meshh>ny)
%         meshh(length(meshh)>ny) = [];
%     end
     dist2GL(1) = dist2GL(2) + calving_front_ice_side;
     dist2GL = dist2GL - dist2GL(1);
      meshh = cumsum(readbin('dely_spacing.bin',[1 ny]));
end
%hydro_cond = [pwd filesep 'sample_CTD_16.mat'];
%hydro_cond = [pwd filesep 'Peter_16.mat'];
%hydro_cond = [pwd filesep 'OD_CTD.mat'];
hydro_cond = [pwd filesep 'openBC.mat'];
if(~isempty(hydro_cond))
    load(hydro_cond);
%    PP = unique(P_nu);
%     for i=1:length(PP)
%         ind = find(PP(i)==P_nu);
%         Salinity(i,1) = mean(Salinity_nu(ind),'omitnan');
%         Temperature(i,1) = mean(Temperature_nu(ind),'omitnan');
%         clear ind;
%     end
%     Pressure = PP;
%     [Pressure,ind] = sort(Pressure,'ascend');
%     Salinity = Salinity(ind);
%     Temperature = Temperature(ind);
end

v0 = 2e3;
h0 = 800;

hfacMin = 0.2;

dlat = 2*0.125/16/16/5; dy=dlat;
dlon = 0.125/1; dx=dlon;

eos = 'jmd95z';
prec = 'real*8';

zC=-delz*([1:nr]-0.5);
zF=-delz*[0:nr];
%size(latc)

% Gravity
gravity= 9.81;
rhoConst= 1030;
rhoIce = 917;

% Nominal depth of model (meters)
H = -500;		%water depth in the ice shelf cavity
Hmin = H + delz;		% deepest point of cavern		
Hmax = H + delz + 150;		% shallowest point of cavern
%jEnd = ny*3/4;		 % where ice-shelf ends
jEnd = ny;
j2=jEnd+1;

if(opt_topo==0)
    long = 0;
    %long=51.4+[0:nx-1]*dlon;
    lonc = long+dlon/2;
    latg = -75.5+[0:ny-1]*dlat;
    %latg=64;
    latc = latg+dlat/2;
    dHdy = (Hmax-Hmin)/dlat/(jEnd-2); %Slope of ice shelf

    bathy = ones(nx,ny)*H;	%For flat bathymetry: bathy = ones(nx,ny)*H;
    bathy(:,1) = 0;

    hIce=Hmin+dHdy*[-1:ny-2]*dlat;
    hIce(1)=0; hIce(j2:ny)=0;
else
    long = 0;
    %long=51.4+[0:nx-1]*dlon;
    lonc = long+dlon/2;
    %latg = -75.5+[0:ny-1]*dlat;
    latg = -75.5+(meshh);
    %latg=64;
    dlat = meshh/deg2km;
    latc = latg+dlat/2;
    
    latc_km = (latc - latc(1)) * deg2km;
    bathy = reshape(interp1(dist2GL,bed_smooth,latc_km),[nx ny]);
    bathy(:,1) = 0;
    ind = min(find(isnan(bathy)));
    if(~isempty(ind))
        bathy(ind:end) = bathy(ind-1);
    end
    
    hIce = reshape(interp1(dist2GL,base_smooth,latc_km),[nx ny]);
    j2 = min(find(isnan(hIce)));
    hIce(1) = 0;
    hIce(j2:ny) = 0;
    H = min(bathy);
    delz = abs(H)/nr;
end
var=([1:ny]-2)/(jEnd-2);
%var(2), var(jEnd)
dMdt_fy=-cos(pi*var); 
%dMdt_fy=0;
dMdt_fy(1)=0; dMdt_fy(j2:ny)=0;

figure(1);clf;
subplot(211)
 yax=[1:ny];
 plot(yax,hIce,'r-');
hold on;
 plot(yax,bathy,'b-');
hold off;
 grid
title('ice-shelf & model depth [m]')

subplot(212)
 plot(yax,dMdt_fy,'b-');
 grid
title('shape of ice-mass input')

namF='bathy_flat.bin';
if kwr > 2,
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,bathy,prec);fclose(fid);
 writebin(namF,bathy);
 fprintf(' done\n');
end

namF='shelficeTopo.Lin.bin';
if kwr > 2,
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,hIce,prec);fclose(fid);
 writebin(namF,hIce);
 fprintf(' done\n');
end

regMsk=ones(ny,1);
regMsk(1)=0; regMsk(j2:ny)=2;
namF='under_Ice_mask.bin';
if kwr > 2,
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,regMsk,prec);fclose(fid);
 writebin(namF,regMsk);
 fprintf(' done\n');
end

%- rate of change due to ice-stream dynamics
%rateDyn=rhoConst*0.1/86400; sfx='r01';
%rateDyn=rhoConst*0.1/3600;  sfx='r02';
 rateDyn=rhoConst*0/3600;  sfx='r02';

dMdt=rateDyn*dMdt_fy;
namF=sprintf('%s.%s.%s','shelfice_dMdt',sfx,'bin');
if kwr > 0,
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,dMdt,prec);fclose(fid);
 writebin(namF,dMdt);
 fprintf(' done\n');
end

dz = delz*ones(1,nr);
zgp1 = [0,cumsum(dz)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);
dz = diff(zgp1);
fprintf('  delZ = %d * %7.6g\n',nr,dz(1))

if(opt_ini==0)
    T_sfc = -0.2;
    T_bot = 0.2;
    del_T = (T_bot - T_sfc)/(59*delz);
    for k = 1:nr;
        tref(k) = T_sfc + del_T*((k-20)*delz);
        tref(k)= max(T_sfc,min(tref(k),T_bot));
    end

    S_sfc = 32;
    S_bot = 34.7;
    del_S = (S_bot - S_sfc)/(59*delz);
    for k = 1:nr;
        sref(k) = S_sfc + del_S*((k-20)*delz);
        sref(k)= max(S_sfc,min(sref(k),S_bot));
    end
else
    % BY VISUAL INSPECTION OF PROFILE
    %Salinity(1:5867) = Salinity(5867);
    %Temperature(1:5867) = Temperature(5867);
    tref = interp1(abs(Pressure),Temperature,abs(zc));
    sref = interp1(abs(Pressure),Salinity,abs(zc));
end
pEOS=-rhoConst*gravity*zC; % in Pa
pEOS=pEOS*1.e-4; % in dBar
rhoAn=densjmd95(sref,tref,pEOS);
rhoAn=rhoAn-rhoConst;

pF=-rhoConst*gravity*zF*1.e-4; % in dBar
rhoUp=densjmd95(sref,tref,pF(2:end));
rhoDw=densjmd95(sref,tref,pF(1:nr));
dRho=rhoUp(1:nr-1)-rhoDw(2:nr);
NSq=-gravity*dRho/delz/rhoConst;

mnV=min(NSq); MxV=max(NSq); Avr=mean(NSq);
fprintf(' Stratif (N^2): min= %e , max= %e , Avr= %e\n',mnV,MxV,Avr);

zax=[1:nr];
% 
% v1=2.5e-2;
% var=1+nr-2*zax; var=var/(nr-1);
% vobc=v1*var;

%% New add to the code (Vobc written so that the flux through the boundary
% is made zero)
% vobc = -0.1/abs(Hmax - Hmin) * (zc + H - (2 * H - Hmax - Hmin)/3);
% vobc(-zc>Hmax) = 0;
% var = 0.1;
% vobc = zeros(size(zc));
% ind = find(zc<abs(H)&zc>abs(H+1));
% vobc(ind) = -var;
% clear ind;
% ind = find(zc<=abs(H+1)&zc>abs(H+3));
% vobc(ind) = (var * abs(zc(ind) - abs(H+1)) - var * abs(zc(ind) - abs(H+3)))/2;
% clear ind;
% ind = find(zc<=abs(H+3)&zc>=Hmax);
% vobc(ind) = var;
vobc = zeros(1,nr);
vobc(1:max(find(zc<abs(hIce(end)+bathy(end))/2))-1) = 0.23/abs(abs(hIce(end)+bathy(end))/2 - abs(bathy(end)));
vobc(1:max(find(zc<abs(hIce(end))))) = 0;
vobc(max(find(zc<abs(hIce(end)+bathy(end))/2))+2:end) = -0.23/abs(abs(hIce(end)+bathy(end))/2 - abs(bathy(end)));
%vobc(end) = -0.1;

figure(2);clf;
 subplot(131);
 plot(tref,-zax,'r-');
 grid
 title('tRef')
 subplot(132);
 plot(sref,-zax,'b-');
 grid
 title('sRef')
 subplot(133);
%plot(vobc,-zax,'b-');
 plot(rhoAn,-zax,'b-');
%plot(NSq*1.e+6,0.5-[2:nr],'r-');
 grid
 title('rhoAn')

if kwr > 2,
 namF='temp_obc.bin';
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,tref,prec);fclose(fid);
 writebin(namF,tref);
 fprintf(' done\n');
 namF='salt_obc.bin';
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,sref,prec);fclose(fid);
 writebin(namF,sref);
 fprintf(' done\n');
 namF='vVel_obc.bin';
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,vobc,prec);fclose(fid);
 writebin(namF,vobc);
 fprintf(' done\n');
end

if kwr > 2,
 var=ones(ny,1)*tref; %size(var)
 namF='temp_ini.bin';
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,var,prec);fclose(fid);
 writebin(namF,var);
 fprintf(' done\n');
%-
 var=ones(ny,1)*sref; %size(var)
 namF='salt_ini.bin';
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,var,prec);fclose(fid);
 writebin(namF,var);
 fprintf(' done\n');
end

rhoAvr=rhoConst-1.345;
fprintf(' convert Ice topo using rhoAve = %10.6f\n',rhoAvr);
mIce0=-rhoAvr*hIce;

namF='shelficeMass.Lin.bin';
if kwr > 1,
 fprintf(' writing file: %s ...',namF);
 %fid=fopen(namF,'w','b'); fwrite(fid,mIce0,prec);fclose(fid);
 writebin(namF,mIce0);
 fprintf(' done\n');
end

return

