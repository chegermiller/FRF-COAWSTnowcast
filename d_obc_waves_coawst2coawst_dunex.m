clear; close all;

t1 = datenum(2021,6,3);
t2 = datenum(2021,6,4);

% Define filenames.
mdl_fname = 'https://geoport.usgs.esipfed.org/thredds/dodsC/vortexfs1/usgs/Projects/COAWST_forecast_refined2/obx_his.ncml';
grd_fname  = 'C:\Users\chegermiller\Desktop\Research\Duck_forecast\models\grids\roms_duck_102019bathy_dunex_small.nc';
init_fname = ['C:\Users\chegermiller\Desktop\Research\Duck_forecast\models\boundary\roms_duck_102019bathy_dunex_small_init_' datestr(t1,'yyyymmdd') '.nc'];
bry_fname  = ['C:\Users\chegermiller\Desktop\Research\Duck_forecast\models\boundary\roms_duck_102019bathy_dunex_small_bry_' datestr(t1,'yyyymmdd') '.nc'];

f_obc_coawst2coawst_dunex(mdl_fname,grd_fname,init_fname,bry_fname,t1,t2)
% f_waves_coawst2coawst_dunex(mdl_fname,t1,t2) % not needed if not running
% ROMS/SWAN configuation
f_FRF2swan_spc2d_dunex(t1,t2)