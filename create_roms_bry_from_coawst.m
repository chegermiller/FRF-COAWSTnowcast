function create_roms_bry_from_coawst(grid_file,bry_file,bry_time,...
    theta_s,theta_b,Tcline,Vtransform,Vstretching,N,...
    zeta_north,ubar_north,vbar_north,...
    ubar_south,vbar_south,zeta_south,...
    zeta_east,ubar_east,vbar_east,...
    zeta_west,ubar_west,vbar_west);

% Create a netcdf file that contains baoundary data for ROMS
% zeta, ubar and vbar.

h=ncread(grid_file,'h');
hmin=min(h(:));
hc=min([hmin,Tcline]);
[LP,MP]=size(h);
L  = LP-1;
M  = MP-1;
xi_psi  = L;
xi_rho  = LP;
xi_u    = L;
xi_v    = LP;
eta_psi = M;
eta_rho = MP;
eta_u   = MP;
eta_v   = M;
%
% These are just copied from above, and then we call get_roms_grid.
%
Sinp.N           =N;            %number of vertical levels
Sinp.Vtransform  =Vtransform;   %vertical transformation equation
Sinp.Vstretching =Vstretching;  %vertical stretching function
Sinp.theta_s     =theta_s;      %surface control parameter
Sinp.theta_b     =theta_b;      %bottom  control parameter
Sinp.Tcline      =Tcline;       %surface/bottom stretching width
Sinp.hc          =hc;           %stretching width used in ROMS
%
Gout=get_roms_grid(grid_file,Sinp);

%create boundary file
create_roms_netcdf_bndry_mwUL(bry_file,Gout,length(bry_time))

%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

ncwrite(bry_file,'zeta_time',bry_time);
ncwrite(bry_file,'v2d_time',bry_time);
ncwrite(bry_file,'v3d_time',bry_time);
ncwrite(bry_file,'salt_time',bry_time);
ncwrite(bry_file,'temp_time',bry_time);

ncwrite(bry_file,'zeta_north',zeta_north);
ncwrite(bry_file,'zeta_south',zeta_south);
ncwrite(bry_file,'zeta_east',zeta_east);
ncwrite(bry_file,'zeta_west',zeta_west);

ncwrite(bry_file,'ubar_north',ubar_north);
ncwrite(bry_file,'ubar_south',ubar_south);
ncwrite(bry_file,'ubar_east',ubar_east);
ncwrite(bry_file,'ubar_west',ubar_west);

ncwrite(bry_file,'vbar_north',vbar_north);
ncwrite(bry_file,'vbar_south',vbar_south);
ncwrite(bry_file,'vbar_east',vbar_east);
ncwrite(bry_file,'vbar_west',vbar_west);

% ncwrite(clm_file,'u',u);
% ncwrite(clm_file,'v',v);
% ncwrite(clm_file,'temp',temp);
% ncwrite(clm_file,'salt',salt);


%close file
disp(['created ', bry_file])


