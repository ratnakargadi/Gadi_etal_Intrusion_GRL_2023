function file_comp(char,timestamp1,timestamp2,nlevels,minvar,maxvar,index,ylims,instname)
%This code will compare the variable in the filenames and plot them as
% a contour plot and finds the difference between the two
addpath('/nfspool-0/home/rgadi/MITgcm/MITgcm-master/utils/matlab');

var1 = squeeze(rdmds(char,timestamp1));
var2 = squeeze(rdmds(char,timestamp2));
Data1 = squeeze(rdmds(instname,timestamp1));
Data2 = squeeze(rdmds(instname,timestamp2));
hFac1 = squeeze(Data1(:,:,end));
hFac2 = squeeze(Data2(:,:,end));
var1(hFac1==0) = NaN;
var2(hFac2==0) = NaN;
% axis_lims = [min(min(var1(:)),min(var2(:))) max(max(var1(:)),max(var2(:)))];
axis_lims = [min(var1(:)) max(var1(:))];
var1 = var1';
var2 = var2';

Y = squeeze(rdmds('YC'));
if(length(unique(Y))==1)
    Y = squeeze(rdmds('XC'));
end

Z = -cumsum(squeeze(rdmds('DRC')));
if(length(Z)~=size(var1,2))
    Z = -cumsum(squeeze(rdmds('DRF')));
end

[Ymesh,Zmesh] = meshgrid(Y,Z);

if(~isempty(minvar))
    [rowmin,colmin] = find(var2<minvar);
else
    rowmin = [];
end

if(~isempty(maxvar))
    [rowmax,colmax] = find(var2>maxvar);
else 
    rowmax = [];
end

figure('units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(2,2,1);
contourf(ax1,Ymesh,Zmesh,var1,nlevels,'LineColor','none');colorbar;
caxis(ax1,axis_lims);
if(~isempty(ylims))
    ylim([ylims]);
end
xlabel(ax1,'Longitude (\circ)');
ylabel(ax1,'Depth (m)');
title(ax1,sprintf('Timestep = %d (a)',timestamp1));

ax2 = subplot(2,2,2);
contourf(ax2,Ymesh,Zmesh,var2,nlevels,'LineColor','none');colorbar;
caxis(ax2,axis_lims);
if(~isempty(ylims))
    ylim([ylims]);
end
xlabel(ax2,'Longitude (\circ)');
ylabel(ax2,'Depth (m)');
title(ax2,sprintf('Timestep = %d (b)',timestamp2));

ax3 = subplot(2,2,3);
vardiff = var2 - var1;
contourf(ax3,Ymesh,Zmesh,vardiff,nlevels,'LineColor','none');colorbar;
caxis(ax3,[-1 1]*max(abs(vardiff(:))))
if(~isempty(ylims))
    ylim([ylims]);
end
if(~isempty(rowmin))
    hold (ax3,'on');
    ymin = Ymesh(min(rowmin),min(colmin));
    zmin = Zmesh(min(rowmin),min(colmin));
    ymax = Ymesh(max(rowmin),max(colmin));
    zmax = Zmesh(max(rowmin),max(colmin));
    rectangle(ax3,'Position',[ymin zmin abs(ymax - ymin) abs(zmax - zmin)],....
        'EdgeColor','k');
end
if(~isempty(rowmax))
    hold (ax3,'on');
    ymin = Ymesh(min(rowmax),min(colmax));
    zmin = Zmesh(min(rowmax),min(colmax));
    ymax = Ymesh(max(rowmax),max(colmax));
    zmax = Zmesh(max(rowmax),max(colmax));
    rectangle(ax3,'Position',[ymin zmin abs(ymax - ymin) abs(zmax - zmin)],....
        'EdgeColor','k');
    hold (ax3,'off');
end
xlim(ax3,[min(Ymesh(:)) max(Ymesh(:))]);
xlabel(ax3,'Longitude (\circ)');
ylabel(ax3,'Depth (m)');
title(ax3,'(b) - (a)');

ax4 = subplot(2,2,4);
if(~isempty(minvar))
    ind = find(var2<minvar);
    var2(ind) = NaN;
end

if(~isempty(maxvar))
    ind = find(var2>maxvar);
    var2(ind) = NaN;
end
vardiff = var2 - var1;
contourf(ax4,Ymesh,Zmesh,vardiff,nlevels,'LineColor','none');colorbar;
caxis(ax4,[-1 1]*max(abs(vardiff(:))))
if(~isempty(ylims))
    ylim([ylims]);
end
xlabel(ax4,'Longitude (\circ)');
ylabel(ax4,'Depth (m)');
title(ax4,'(b) - (a)');

switch char
    case 'S'
        str = 'Salinity';
    case 'T'
        str = 'Temperature';
    case 'U'
        str = 'Zonal velocity';
    case 'V'
        str = 'Meridonal velocity';
    case 'W'
        str = 'vertical velocity';
    case 'Eta'
        str = 'Ice shelf interface deviation';
    otherwise
        str = 'unknown';
end
      
suptitle(sprintf('%s',str));
str2 = [str '_' num2str(timestamp1,'%d') '_' num2str('timestamp2','%d') '_' num2str(index,'%d')];
print(gcf,'-dpng','-r300',fullfile(pwd,str2));
end

