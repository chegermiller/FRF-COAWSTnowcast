function f_FRF2swan_spc2d_dunex(t1,t2)
% go from FRF 8 m directional wave array .nc to SWAN .spc2d forcing file
%
% https://chlthredds.erdc.dren.mil/thredds/catalog/frf/oceanography/waves/8m-array/catalog.html
%
% 4/2019 CAHegermiller

url = ['https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/waves/8m-array/2021/FRF-ocean_waves_8m-array_' datestr(t1,'yyyymm') '.nc'];
txt_fname = ['C:\Users\chegermiller\Desktop\Research\Duck_forecast\models\boundary\8m_array_' datestr(t1,'yyyymmdd') '.spc2d'];

%% load spc2d nc files
time = ncread(url,'time')/86400+datenum(1970,1,1);
tstart = find(time == t1);
tend = find(time == t2);
time = time(tstart:tend);

lon = ncread(url,'longitude');
lat = ncread(url,'latitude');

f = ncread(url,'waveFrequency');
d = ncread(url,'waveDirectionBins');

E = ncread(url,'directionalWaveEnergyDensity',[1 1 tstart],[Inf Inf tend-tstart+1]);

E(isnan(E)) = -0.9900e2;
E(E>=9.9999e-7 & E<=1.0000e-6)=0;
E = E * 1025*9.81; 

nf = length(f);
nd = length(d);

% datenum to swan date format
dn = datestr(time,'yyyymmdd.HHMMSS');

% header
fid = fopen(txt_fname,'w');
fprintf(fid,'%s\n','SWAN   1                                Swan standard spectral file, version');
fprintf(fid,'%s\n','$   Data produced by SWAN version 41.20');
fprintf(fid,'%s\n','$   Project: DUNEX    ;  run number:');
fprintf(fid,'%s\n','TIME                                    time-dependent data');
fprintf(fid,'%s\n','     1                                  time coding option');
fprintf(fid,'%s\n','LONLAT                                  locations in spherical coordinates');
fprintf(fid,'%s\n','     1                                  number of locations');
fprintf(fid,'   %.6f   %.6f\n',lon,lat);
fprintf(fid,'%s\n','AFREQ                                   absolute frequencies in Hz');
fprintf(fid,'    %d                                  number of frequencies\n',nf);
fprintf(fid,'    %.4f\n',f');
fprintf(fid,'%s\n','NDIR                                    spectral nautical directions in degr');
fprintf(fid,'    %d                                  number of directions\n',nd);
fprintf(fid,'  %.4f\n',d);
fprintf(fid,'%s\n','QUANT');
fprintf(fid,'%s\n','     1                                  number of quantities in table');
fprintf(fid,'%s\n','EnDens                                  energy densities in J/m2/Hz/degr');
fprintf(fid,'%s\n','J/m2/Hz/degr                            unit');
fprintf(fid,'%s\n','   -0.9900E+02                          exception value');

for ii = 1:length(dn)
    spec = sq(E(:,:,ii));
    spec = spec';
    mins = min(spec(spec~=0));
    mult = 1/mins;
    spec = round(spec*mult);
  
    fprintf(fid,'%s                         date and time\n',dn(ii,:));
    fprintf(fid,'%s\n','FACTOR');
    fprintf(fid,'    %.8E\n',mins);
    for ff = 1:nf
        fprintf(fid,[repmat('    %d',1,nd) '\n'],spec(ff,:));
    end
end

fclose(fid);
end