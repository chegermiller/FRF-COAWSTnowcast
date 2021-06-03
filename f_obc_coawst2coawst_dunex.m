function f_obc_coawst2coawst_dunex(mdl_fname,grd_fname,init_fname,bry_fname,t1,t2)
% This routine :
%  - creates boundary and initial condition files for ROMS:
%    coawst_bdy.nc ; coawst_ini.nc
%    on a user-defined grid for a user-defined date.
%
% bry file includes zeta, ubar, vbar
% ini file includes zeta, u, v, ubar, vbar, temp, salt
%
% This is currently set up to use opendap and nctoolbox.
%
% written by Maitane Olabarrieta 01/14/2018
% CAH edits

% Enter grid vertical coordinate parameters
% These need to be consistent with the ROMS setup.
theta_s = 0.0;
theta_b = 0.0;
Tcline  = 20;
N       = 5;
Vtransform  = 2;
Vstretching = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('getting roms grid dimensions ...');

Sinp.N           = N;            % number of vertical levels
Sinp.Vtransform  = Vtransform;   % vertical transformation equation
Sinp.Vstretching = Vstretching;  % vertical stretching function
Sinp.theta_s     = theta_s;      % surface control parameter
Sinp.theta_b     = theta_b;      % bottom  control parameter
Sinp.Tcline      = Tcline;       % surface/bottom stretching width

if Vtransform == 1
    h = ncread(grd_fname,'h');
    hmin = min(h(:));
    hc = min(max(hmin,0),Tcline);
elseif Vtransform == 2
    h = ncread(grd_fname,'h');
    hmin = max(0.1,min(h(:)));
    hc = Tcline;
end

Sinp.hc = hc;                    % stretching width used in ROMS
gn = get_roms_grid(grd_fname,Sinp);
[nxr,nyr] = size(gn.lon_rho);
[nxu,nyu] = size(gn.lon_u);
[nxv,nyv] = size(gn.lon_v);

%% Read data

setup_nctoolbox;
nc = ncgeodataset(mdl_fname);
time = nj_time(nc,'zeta');

tstart = find(time == t1);
tend = find(time == t2);

time = time(tstart:tend);
T = length(time);

% GET HORIZONTAL GRID
g = nc{'zeta'}(1,:,:).grid;
lon_rho = g.lon;
lat_rho = g.lat;

% READ GRID PARAMETERS
masku = nc{'mask_u'}(:);
maskv = nc{'mask_v'}(:);
maskr = nc{'mask_rho'}(:);
angler = nc{'angle'}(:);
coawst_h = nc{'h'}(:);
coawst_Vtransform = nc{'Vtransform'}(:);
coawst_Vstretching = nc{'Vstretching'}(:);
coawst_hc = nc{'hc'}(:);
coawst_theta_s = nc{'theta_s'}(:);
coawst_theta_b = nc{'theta_b'}(:);
coawst_N = length(nc{'s_rho'}(:));

lon_rho_z = repmat(lon_rho,1,1,coawst_N);
lat_rho_z = repmat(lat_rho,1,1,coawst_N);

% READ VARIABLES
count = 1;
for mm = tstart:tend
    zeta_coawst(:,:,count)   = double(sq(nc{'zeta'}(mm,:,:)));
    ubar_coawst(:,:,count)   = double(sq(nc{'ubar'}(mm,:,:)));
    vbar_coawst(:,:,count)   = double(sq(nc{'vbar'}(mm,:,:)));
    wmasku(:,:,count) = nc{'wetdry_mask_u'}(mm,:,:);
    wmaskv(:,:,count) = nc{'wetdry_mask_v'}(mm,:,:);
    wmaskr(:,:,count) = nc{'wetdry_mask_rho'}(mm,:,:);
end

mm = tstart;
u_coawst    = double(sq(nc{'u'}(mm,:,:,:)));
v_coawst    = double(sq(nc{'v'}(mm,:,:,:)));
temp_coawst = double(sq(nc{'temp'}(mm,:,:,:)));
salt_coawst = double(sq(nc{'salt'}(mm,:,:,:)));

save coawst_vars.mat zeta_coawst ubar_coawst vbar_coawst u_coawst v_coawst temp_coawst salt_coawst -v7.3

close(nc);

%% Initialize the variables of interest

u = zeros(nxu,nyu,N);
v = zeros(nxv,nyv,N);
temp = zeros(nxr,nyr,N);
salt = zeros(nxr,nyr,N);
ubar = zeros(nxu,nyu,T);
vbar = zeros(nxv,nyv,T);
zeta = zeros(nxr,nyr,T);

lon_rho_romsz = repmat(gn.lon_rho,1,1,N);
lat_rho_romsz = repmat(gn.lat_rho,1,1,N);

%% Compute depth at u, v points for both grids and get rid of land or masked cells

% hu = rho2u_2d_mw(h);
% hu(~gn.mask_u) = NaN;
% hu(hu<=0) = NaN;
% hv = rho2v_2d_mw(h);
% hv(~gn.mask_v) = NaN;
% hv(hv<=0) = NaN;
%  
% coawst_hu = rho2u_2d_mai(coawst_h);
% coawst_hu(~masku) = NaN;
% coawst_hu(coawst_hu<=0) = NaN;
% coawst_hv = rho2v_2d_mai(coawst_h);
% coawst_hv(~maskv) = NaN;
% coawst_hv(coawst_hv<=0) = NaN;

%% READ 3D & 4D VARIABLES FOR INIT

for mm = 1
    aa = sq(zeta_coawst(:,:,mm));
    
    % Compute vertical elevations of the grid, this is time dependent
    % These depths are negative
    if mod(mm-1,24) == 0
        zr = set_depth(coawst_Vtransform,coawst_Vstretching,coawst_theta_s,coawst_theta_b,coawst_hc,coawst_N, ...
            1,coawst_h,aa);
    end
    
    aa(~maskr) = NaN;
    aa = maplev(aa);
    zz = griddata(lon_rho(:),lat_rho(:),aa(:),gn.lon_rho,gn.lat_rho);
    test = zz;
    S = size(gn.lon_rho);
    sigma = 5;
    test = conv2(test,gaussian2d(S,sigma),'same');
    zz(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3)) = test(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3));
    zeta(:,:,mm) = zz;
    
    au = sq(ubar_coawst(:,:,mm));
    au(~masku) = NaN;
    au = maplev(au);
    % au = au.*coawst_hu;
    aur = u2rho_2d_mai(au);
    
    av = sq(vbar_coawst(:,:,mm));
    av(~maskv) = NaN;
    av = maplev(av);
    % av = av.*coawst_hv;
    avr = v2rho_2d_mai(av);
    
    % Compute northward and eastward velocities, important!
    vel = aur + avr.*sqrt(-1);
    vel = vel .* exp(sqrt(-1)*angler);
    velu = real(vel);
    velv = imag(vel);
    velu1 = griddata(lon_rho(:),lat_rho(:),velu(:),gn.lon_rho,gn.lat_rho);
    velv1 = griddata(lon_rho(:),lat_rho(:),velv(:),gn.lon_rho,gn.lat_rho);
    test = velu1;
    test = conv2(test,gaussian2d(S,sigma),'same');
    velu1(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3)) = test(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3));
    test = velv1;
    test = conv2(test,gaussian2d(S,sigma),'same');
    velv1(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3)) = test(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3));
    
    % Rotate velocities to ROMS grid, important!
    ubar1 = velu1.*cos(gn.angle)+velv1.*sin(gn.angle);
    vbar1 = velv1.*cos(gn.angle)-velu1.*sin(gn.angle);
    ubar(:,:,mm)=rho2u_2d_mw(ubar1);  % defined at u points
    vbar(:,:,mm)=rho2v_2d_mw(vbar1); 
    
%     ubar1 = rho2u_2d_mw(ubar1);
%     vbar1 = rho2v_2d_mw(vbar1);
%     ubar1 = ubar1./hu;
%     vbar1 = vbar1./hv;
%     ubar(:,:,mm) = maplev(ubar1);
%     vbar(:,:,mm) = maplev(vbar1);
    
    clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
    
    for zz = 1:coawst_N
        aa = sq(temp_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        temp2(:,:,zz) = aa;
        clear aa
        
        aa = sq(salt_coawst(zz,:,:));
        aa(~maskr) = NaN;
        aa = maplev(aa);
        salt2(:,:,zz) = aa;
        clear aa
        
        aa = sq(u_coawst(zz,:,:));
        aa(~masku) = NaN;
        aa = maplev(aa);
        u2(:,:,zz) = aa;
        clear aa
        
        aa = sq(v_coawst(zz,:,:));
        aa(~maskv) = NaN;
        aa = maplev(aa);
        v2(:,:,zz) = aa;
        clear aa
    end
    
    for zz = 1:coawst_N
        zr(:,:,zz) = maplev(sq(zr(:,:,zz)));
    end
    
    temp = griddata(lon_rho_z,lat_rho_z,zr,temp2,lon_rho_romsz,lat_rho_romsz,gn.z_r);
    salt = griddata(lon_rho_z,lat_rho_z,zr,salt2,lon_rho_romsz,lat_rho_romsz,gn.z_r);
    
    for zz = 1:coawst_N
        ur(:,:,zz) = u2rho_2d_mai(u2(:,:,zz));
        vr(:,:,zz) = v2rho_2d_mai(v2(:,:,zz));
    end
    
    % Compute Northward and Eastward velocities
    for zz = 1:coawst_N
        vel = sq(ur(:,:,zz))+sq(vr(:,:,zz)).*sqrt(-1);
        vel = vel.* exp(sqrt(-1) * angler);
        ur(:,:,zz) = real(vel);
        vr(:,:,zz) = imag(vel);
    end
    
    ur2 = griddata(lon_rho_z,lat_rho_z,zr,ur,lon_rho_romsz,lat_rho_romsz,gn.z_r);
    vr2 = griddata(lon_rho_z,lat_rho_z,zr,vr,lon_rho_romsz,lat_rho_romsz,gn.z_r);
    
    % Rotate velocities to ROMS grid, important!
    for zz = 1:N
        ur2(:,:,zz) = sq(ur2(:,:,zz)).*cos(gn.angle)+sq(vr2(:,:,zz)).*sin(gn.angle);
        vr2(:,:,zz) = sq(vr2(:,:,zz)).*cos(gn.angle)-sq(ur2(:,:,zz)).*sin(gn.angle);
        u(:,:,zz) = rho2u_2d_mw(ur2(:,:,zz));  % defined at u points
        v(:,:,zz) = rho2v_2d_mw(vr2(:,:,zz));  % defined at v points
    end
    
    % Remove possible NaN values
    for zz = 1:N
        aa = sq(temp(:,:,zz));
        aa = maplev(aa);
        test = aa;
        S = size(gn.lon_rho);
        test = conv2(test,gaussian2d(S,sigma),'same');
        aa(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3)) = test(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3));
        temp(:,:,zz) = aa;
        clear aa
        
        aa = sq(salt(:,:,zz));
        aa = maplev(aa);
        test = aa;
        test = conv2(test,gaussian2d(S,sigma),'same');
        aa(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3)) = test(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3));
        salt(:,:,zz) = aa;
        clear aa
        
        aa = sq(u(:,:,zz));
        aa = maplev(aa);
        test = aa;
        S = size(gn.lon_u);
        test = conv2(test,gaussian2d(S,sigma),'same');
        aa(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3)) = test(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3));
        u(:,:,zz) = aa;
        clear aa
        
        aa = sq(v(:,:,zz));
        aa = maplev(aa);
        test = aa;
        S = size(gn.lon_v);
        test = conv2(test,gaussian2d(S,sigma),'same');
        aa(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3)) = test(sigma*3:end-(sigma*3),sigma*3:end-(sigma*3));
        v(:,:,zz) = aa;
        clear aa
    end
    
    %% CREATE INITIAL CONDITION
    init_time = time(1)-datenum(1858,11,17);
    create_roms_init_from_coawst(grd_fname,init_fname,init_time,...
        Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,...
        Sinp.N,u,v,sq(ubar(:,:,mm)),sq(vbar(:,:,mm)),...
        temp,salt,sq(zeta(:,:,mm)));
    % little cheat
    ncwrite(init_fname,'hc',2e1);
end
   
%% READ 3D VARIABLES FOR BRY

for mm = 2:T
    
    aa = sq(zeta_coawst(:,:,mm));
    
    % Compute vertical elevations of the grid, this is time dependent
    % These depths are negative
    
%     if mod(mm-1,24) == 0
%         zr = set_depth(coawst_Vtransform,coawst_Vstretching,coawst_theta_s,coawst_theta_b,coawst_hc,coawst_N, ...
%             1,coawst_h,aa);
%     end
    
    aa(~sq(wmaskr(:,:,mm))) = NaN;
    % aa(~maskr) = NaN;
    %aa = maplev(aa);
    zz = griddata(lon_rho(:),lat_rho(:),aa(:),gn.lon_rho,gn.lat_rho);
    zz = maplev(zz);
    zeta(:,:,mm) = zz;
    
    au = sq(ubar_coawst(:,:,mm));
    au(~sq(wmasku(:,:,mm))) = NaN;
    %au(~masku) = NaN;
    %au = maplev(au);
    %au = au.*coawst_hu;
    aur = u2rho_2d_mai(au);
    
    av = sq(vbar_coawst(:,:,mm));
    av(~sq(wmaskv(:,:,mm))) = NaN;
    %av(~maskv) = NaN;
    %av = maplev(av);
    %av = av.*coawst_hv;
    avr = v2rho_2d_mai(av);
    
    % Compute northward and eastward velocities, important!
    vel = aur + avr.*sqrt(-1);
    vel = vel .* exp(sqrt(-1)*angler);
    velu = real(vel);
    velv = imag(vel);
    velu1 = griddata(lon_rho(:),lat_rho(:),velu(:),gn.lon_rho,gn.lat_rho);
    velv1 = griddata(lon_rho(:),lat_rho(:),velv(:),gn.lon_rho,gn.lat_rho);
    
    % Rotate velocities to ROMS grid, important!
    ubar1 = velu1.*cos(gn.angle)+velv1.*sin(gn.angle);
    vbar1 = velv1.*cos(gn.angle)-velu1.*sin(gn.angle);
    ubar1 = maplev(ubar1);
    vbar1 = maplev(vbar1);
    ubar(:,:,mm)=rho2u_2d_mw(ubar1);  % defined at u points
    vbar(:,:,mm)=rho2v_2d_mw(vbar1);
    
%     ubar1 = rho2u_2d_mw(ubar1);
%     vbar1 = rho2v_2d_mw(vbar1);
%     ubar1 = ubar1./hu;
%     vbar1 = vbar1./hv;
    
    clear vel au av aur avr velu velv velu1 velv1 ubar1 vbar1
end

%% CREATE BOUNDARY CONDITION

%time = time(ind:end);
time = time-datenum(1858,11,17);

zeta_north = sq(zeta(:,end,:));
ubar_north = sq(ubar(:,end,:));
vbar_north = sq(vbar(:,end,:));

zeta_south = sq(zeta(:,1,:));
ubar_south = sq(ubar(:,1,:));
vbar_south = sq(vbar(:,1,:));

zeta_east = sq(zeta(end,:,:));
ubar_east = sq(ubar(end,:,:));
vbar_east = sq(vbar(end,:,:));

zeta_west = sq(zeta(1,:,:));
ubar_west = sq(ubar(1,:,:));
vbar_west = sq(vbar(1,:,:));

create_roms_bry_from_coawst(grd_fname,bry_fname,time,...
Sinp.theta_s,Sinp.theta_b,Sinp.Tcline,Sinp.Vtransform,Sinp.Vstretching,Sinp.N,...
zeta_north,ubar_north,vbar_north,...
ubar_south,vbar_south,zeta_south,...
zeta_east,ubar_east,vbar_east,...
zeta_west,ubar_west,vbar_west);
end