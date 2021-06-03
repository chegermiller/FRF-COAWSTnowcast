function f_waves_coawst2coawst_dunex(url,t1,t2)
% This routine :
%  - creates TPAR files for SWAN boundary forcing at pre-defined boundary
%    points on a domain
%
% This is currently set up to use opendap and nctoolbox.
%
% CAH

setup_nctoolbox;

load('C:\Users\chegermiller\Desktop\Research\Duck_forecast\models\grids\specpts_duck_dunex.mat');

nc = ncgeodataset(url);
time = nj_time(nc,'zeta');

tstart = find(time == t1);
tend = find(time == t2);

xg = ncread(url,'lon_rho');
yg = ncread(url,'lat_rho');

time = time(tstart:tend);
T = length(time);

Hwave = nc{'Hwave'}(tstart:tend,:,:);
Hwave = permute(Hwave,[3,2,1]);
Pwave = nc{'Pwave_top'}(tstart:tend,:,:);
Pwave = permute(Pwave,[3,2,1]);
Dwave = nc{'Dwave'}(tstart:tend,:,:);
Dwave = permute(Dwave,[3,2,1]);

% this is a hack to use 8m array directional spread since we don't output
url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/waves/8m-array/2021/FRF-ocean_waves_8m-array_202104.nc';
time2 = ncread(url,'time')/86400+datenum(1970,1,1);

tstart = find(time2 == t1);
tend = find(time2 == t2);

spread = ncread(url,'directionalPeakSpread');
spread = interp1(time2,spread,time);

for ss = 1:length(specpts)
    gx = specpts(ss,1);
    gy = specpts(ss,2);
    
    [idx] = nearxy(xg(:),yg(:),gx,gy);
    [ii,jj] = ind2ij(xg,idx);
    
    daplon = xg(ii-1:ii+1,jj-1:jj+1);
    daplat = yg(ii-1:ii+1,jj-1:jj+1);
    
    TPAR(:,1) = str2num(datestr(time,'yyyymmdd.HHMM'));
    
    temp = Hwave(ii-1:ii+1,jj-1:jj+1,:);
    temp(temp>1000) = 0;
    for tt = 1:length(time)
        test = sq(temp(:,:,tt));
        F = scatteredInterpolant(daplon(:),daplat(:),double(test(:)),'linear','none');
        cff = F(gx,gy);
        cff(isnan(cff)) = 0;
        TPAR(tt,2) = cff;
    end
    
    temp = Pwave(ii-1:ii+1,jj-1:jj+1,:);
    temp(temp>1000) = 0;
    for tt = 1:length(time)
        test = sq(temp(:,:,tt));
        F = scatteredInterpolant(daplon(:),daplat(:),double(test(:)),'linear','none');
        cff = F(gx,gy);
        cff(isnan(cff)) = 0;
        TPAR(tt,3) = cff;
    end
    
    temp = Dwave(ii-1:ii+1,jj-1:jj+1,:);
    temp(temp>1000) = 0;
    for tt = 1:length(time)
        test = sq(temp(:,:,tt));
        F = scatteredInterpolant(daplon(:),daplat(:),double(test(:)),'linear','none');
        cff = F(gx,gy);
        cff(isnan(cff)) = 0;
        TPAR(tt,4) = cff;
    end
    
    TPAR(:,5) = spread(ii);
    
    fname = ['TPAR',num2str(ss),'.txt'];
    fid = fopen(fname,'w');
    fprintf(fid,'TPAR \n');
    fprintf(fid,'%8.4f         %3.2f        %3.2f     %3.f.       %2.f\n',TPAR');
    fclose(fid);
end

%DirIOfname = 'C:\Users\chegermiller\Desktop\Research\Duck_forecast\models\grids\specpts_duck_dunex.mat';
%runpath = '';
%bdry_com(DirIOfname,runpath);
end