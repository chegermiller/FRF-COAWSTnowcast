# FRF-COAWSTnowcast
This repo includes scripts used to generated boundary and initial conditions for a high resolution COAWST (ROMS/InWave) application at the FRF, Duck, NC, USA.

d_obc_waves_coawst2coawst_dunex.m 
Driver script to create open boundary conditions (obc) and wave boundary conditions for a ROMS/InWave application.
This script will call the functions:

f_obc_coawst2coawst_dunex.m
This function extracts variables from the COAWST forecast for a user-defined time period and interpolates them to a user-defined grid. 
An initial condition file (.nc) and an open boundary condition file (.nc) with user-defined names will be generated.
Requires create_roms_init_from_coawst.m and create_roms_bry_from_coawst.m.

f_waves_coawst2coawst_dunex.m
This function extracts bulk waves parameters from the COAWST forecast for a user-defined time period and writes SWAN TPAR files (.txt) at predefined boundary points on a user-defined grid.
Commented out because this is not necessary if not running a ROMS/SWAN application.

f_FRF2swan_spc2d_dunex.m
This function pulls directional wave spectra at the 8 m array from the FRF THREDDS server for a user-defined time period and write a SWAN-style spc2d file (.txt).

grid = roms_duck_102019bathy_dunex_small.nc
example obc = roms_duck_102019bathy_dunex_small_bry_20210509.nc
example init = roms_duck_102019bathy_dunex_small_init_20210509.nc
example spc2d = 8m_array_20210509.spc2d
