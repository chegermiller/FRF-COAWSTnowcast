function create_roms_init_from_coawst(grid_file,init_file,init_time,...
theta_s,theta_b,Tcline,Vtransform,Vstretching,N,u1,v1,ubar1,vbar1,...
temp1,salt1,zeta1)

% Create a netcdf file that contains initialization data for ROMS.
% Initializes temp, salt, u, v, ubar, vbar, and 
% all sediment parameters.
%
% In cppdefs.h you should have 
% #undef ana_initial
% #undef ana_sediment
%
% Users can adapt this file to their own application.
%
% jcw 7-8-2008
% jcw 07Feb2012 updated for netcdf mathworks interface
% jcw 19March2015 update to use get_roms_grid
% jcw 28July2016 updated veg params
% mo 01Jan2018 modified to initialized from coawst operational
%

%!         W-level  RHO-level                                           !
%!                                                                      !
%!            N     _________                                           !
%!                 |         |                                          !
%!                 |    N    |                                          !
%!          N-1    |_________|  S                                       !
%!                 |         |  E                                       !
%!                 |   N-1   |  A                                       !
%!            2    |_________|  W                                       !
%!                 |         |  A                                       !
%!                 |    2    |  T                                       !
%!            1    |_________|  E                                       !
%!                 |         |  R                                       !
%!                 |    1    |                                          !
%!            0    |_________|_____ bathymetry                          !
%!                 |/////////|                                          !
%!                 |    1    |                                          !
%!            1    |_________|  S                                       !
%!                 |         |  E                                       !
%!                 |    2    |  D                                       !
%!            2    |_________|  I                                       !
%!                 |         |  M                                       !
%!                 |  Nbed-1 |  E                                       !
%!        Nbed-1   |_________|  N                                       !
%!                 |         |  T                                       !
%!                 |  Nbed   |                                          !
%!         Nbed    |_________|                                          !
%!                                                                      !

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc some grid stuff here - do not change this.
% You go on to step 6.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
%6) Initialize zeta, salt, temp, u, v, ubar, vbar.
%
% Init values for zeta.
  display('Initializing zeta')
  zeta(1:xi_rho,1:eta_rho,1:length(init_time)) = zeta1;

%
% Init values for u, ubar, v, and vbar.
  display('Initializing u, v, ubar, and vbar')

    u(1:xi_u,1:eta_u,1:N,1:length(init_time)) = u1;
    ubar(1:xi_u,1:eta_u,1:length(init_time)) = ubar1;
    v(1:xi_v,1:eta_v,1:N,1:length(init_time)) = v1;
    vbar(1:xi_v,1:eta_v,1:length(init_time)) = vbar1;
%
% Init values for temp and salt.
  display('Initializing temp and salt')
%
  NAT=2;                            % Number of active tracers. Usually 2 (temp + salt).
  salt(1:xi_rho,1:eta_rho,1:N,1:length(init_time)) = salt1;
  temp(1:xi_rho,1:eta_rho,1:N,1:length(init_time)) = temp1;

% %7) Enter number of mud sediments (NCS) and number of sand sediments (NNS).
% %   These values should be the same as in mod_param.F
%     NCS = 0;   %number of cohesive sed classes
%     NNS = 1;   %number of non-cohesive sed classes
% %
% % calc sed parameters. Do not alter.
% %
%    NST = NCS + NNS;     % total number of sed tracers.
%    NT = NAT+NST;        % total number of tracers.
% 
% %8) Enter number of bed layers
% %   This value should be the same as in mod_param.F
%     Nbed = 1;
% 
% %9) Sediment class properties (in order, mud first then sand).
% %  These values should coincide with your sediment.in file.
%   mud_Srho=ones(1,NCS)*2650;        %kg m-3, NCS values
%   mud_Sd50=[0.06]/1000;             %m,      NCS values
%   mud_Wsed=[1.0]/1000;              %m s-1,  NCS values
%   mud_tau_ce=[0.05];                %N m-2,  NCS values
%   mud_Erate=[5]*1e-5;               %kg m-2 s-1, NCS values
%   sand_Srho=ones(1,NNS)*2650;       %kg m-3, NNS values
%   sand_Sd50=[0.2]/1000;             %m,      NNS values
%   sand_Wsed=[27]/1000;              %m s-1,  NNS values
%   sand_tau_ce=[0.19];               %N m-2,  NNS values
%   sand_Erate=[1]*1e-5;              %kg m-2 s-1, NNS values
% %
% % make some combined arrays.  Do not alter.
% %
%   Srho=  [mud_Srho,sand_Srho];
%   Sd50=  [mud_Sd50,sand_Sd50];
%   Wsed=  [mud_Wsed,sand_Wsed];
%   tau_ce=[mud_tau_ce,sand_tau_ce];
%   Erate= [mud_Erate,sand_Erate];
% 
% %9) Provide initial sediment properties in water column.
%   display('Initializing suspended sediments.')
% %
% % mud.
% %
% for idmud=1:NCS
%   count=['0',num2str(idmud)];
%   count=count(end-1:end);
%   zz=zeros(xi_rho,eta_rho,N);
%   eval(['mud_',count,'(1:xi_rho,1:eta_rho,1:N,1:length(init_time)) = zz;'])               %mud conc in water column
% end
% %
% % sand.
% %
% for isand=1:NNS
%   count=['0',num2str(isand)];
%   count=count(end-1:end);
%   zz=zeros(xi_rho,eta_rho,N);
%   eval(['sand_',count,'(1:xi_rho,1:eta_rho,1:N,1:length(init_time)) = zz;'])               %sand conc in water column
% end
% 
% %10) Provide initial sediment properties in bed.
% %
% % bed properties
%     display('Initializing sediment bed.')
%     bed_thickness(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time)) = 1.0;
%     bed_age(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time))       = init_time(1);
%     bed_porosity(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time))  = 0.5;
%     bed_biodiff(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time))   = 0.0;
% %
% % for mud
% %
% for idsed=1:NCS
%   count=['0',num2str(idsed)];
%   count=count(end-1:end);
%   eval(['mudfrac_',count,'(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time)) = 1/NST;'])      %fraction of each sed class in each bed cell
%   eval(['mudmass_',count,'(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time)) = bed_thickness(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time)).*Srho(idsed).*(1.0-bed_porosity(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time))).*mudfrac_',count,'(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time));'])          %mass of each sed class in each bed cell
% end
% %
% % for sand
% %
% for idsed=1:NNS
%  count=['0',num2str(idsed)];
%  count=count(end-1:end);
%  eval(['sandfrac_',count,'(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time)) = 1/NST;'])      %fraction of each sed class in each bed cell
%  eval(['sandmass_',count,'(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time)) = bed_thickness(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time)).*Srho(idsed).*(1.0-bed_porosity(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time))).*sandfrac_',count,'(1:xi_rho,1:eta_rho,1:Nbed,1:length(init_time));'])          %mass of each sed class in each bed cell
%  eval(['bedload_Usand_',count,'(1:xi_u,1:eta_u,1:length(init_time)) = 0;'])              %bed load
%  eval(['bedload_Vsand_',count,'(1:xi_v,1:eta_v,1:length(init_time)) = 0;'])              %bed load
% end
% 
% %11)
% %
% % set some surface properties
%   display('Initializing sediment surface properties.')
% %
% cff1=1.0;
% cff2=1.0;
% cff3=1.0;
% cff4=1.0;
% for ised=1:NCS
%   count=['0',num2str(ised)];
%   count=count(end-1:end);
%   eval(['cff1=cff1.*mud_Sd50(ised).^squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
%   eval(['cff2=cff2.*mud_Srho(ised).^squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
%   eval(['cff3=cff3.*mud_Wsed(ised).^squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
%   eval(['cff4=cff4.*mud_tau_ce(ised).^squeeze(mudfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
% end
% for ised=1:NNS
%   count=['0',num2str(ised)];
%   count=count(end-1:end);
%   eval(['cff1=cff1.*sand_Sd50(ised).^squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
%   eval(['cff2=cff2.*sand_Srho(ised).^squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
%   eval(['cff3=cff3.*sand_Wsed(ised).^squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
%   eval(['cff4=cff4.*sand_tau_ce(ised).^squeeze(sandfrac_',count,'(1:xi_rho,1:eta_rho,1,1));'])
% end
% grain_diameter=cff1;
% grain_density=cff2;
% settling_vel=cff3;
% erosion_stress=cff4;
% ripple_length=0.10;
% ripple_height=0.01;
% dmix_offset=0.0;
% dmix_slope=0.0;
% dmix_time=0.0;
% 
% %12)
% %
% % set vegetation properties 
% NVEG=1;
% plant_density=zeros(xi_rho,eta_rho,NVEG);
% plant_height=zeros(xi_rho,eta_rho,NVEG);
% plant_diameter=zeros(xi_rho,eta_rho,NVEG);
% plant_thickness=zeros(xi_rho,eta_rho,NVEG);
% marsh_mask=zeros(xi_rho,eta_rho,NVEG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  END of USER INPUT                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create init file
Nbed = 0;
NNS = 0;
NCS = 0;
NVEG = 0;
create_roms_netcdf_init_mw(init_file,Gout,Nbed,NNS,NCS,NVEG)

%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

ncwrite(init_file,'theta_s',Gout.theta_s);
ncwrite(init_file,'theta_b',Gout.theta_b);
ncwrite(init_file,'Tcline',Gout.Tcline);
ncwrite(init_file,'Cs_r',Gout.Cs_r);
ncwrite(init_file,'Cs_w',Gout.Cs_w);
ncwrite(init_file,'sc_w',Gout.s_w);
ncwrite(init_file,'sc_r',Gout.s_rho);
ncwrite(init_file,'hc',Gout.hc);
ncwrite(init_file,'Vtransform',Gout.Vtransform);
ncwrite(init_file,'Vstretching',Gout.Vstretching);
ncwrite(init_file,'spherical',Gout.spherical);

ncwrite(init_file,'ocean_time',init_time);

ncwrite(init_file,'zeta',zeta);
ncwrite(init_file,'ubar',ubar);
ncwrite(init_file,'vbar',vbar);
ncwrite(init_file,'u',u);
ncwrite(init_file,'v',v);
ncwrite(init_file,'temp',temp);
ncwrite(init_file,'salt',salt);

% for mm=1:NCS
%   count=['00',num2str(mm)];
%   count=count(end-1:end);
%   eval(['ncwrite(init_file,''mud_',count,''',mud_',count,');'])            %sed conc in water column
%   eval(['ncwrite(init_file,''mudfrac_',count,''',mudfrac_',count,');'])    %sed frac on bed
%   eval(['ncwrite(init_file,''mudmass_',count,''',mudmass_',count,');'])    %sed mass on bed
% end
% for mm=1:NNS
%   count=['00',num2str(mm)];
%   count=count(end-1:end);
%   eval(['ncwrite(init_file,''sand_',count,''',sand_',count,');'])           %sed conc in water column
%   eval(['ncwrite(init_file,''sandfrac_',count,''',sandfrac_',count,');'])   %sed frac on bed
%   eval(['ncwrite(init_file,''sandmass_',count,''',sandmass_',count,');'])   %sed mass on bed
% %
%   eval(['ncwrite(init_file,''bedload_Usand_',count,''',bedload_Usand_',count,');'])   %bedload
%   eval(['ncwrite(init_file,''bedload_Vsand_',count,''',bedload_Vsand_',count,');'])   %bedload
% end
% 
% ncwrite(init_file,'bed_thickness',bed_thickness);
% ncwrite(init_file,'bed_age',bed_age);
% ncwrite(init_file,'bed_porosity',bed_porosity);
% ncwrite(init_file,'bed_biodiff',bed_biodiff);
% 
% ncwrite(init_file,'grain_diameter',grain_diameter);
% ncwrite(init_file,'grain_density',grain_density);
% ncwrite(init_file,'settling_vel',settling_vel);
% ncwrite(init_file,'erosion_stress',erosion_stress);
% ncwrite(init_file,'ripple_height',ripple_height);
% ncwrite(init_file,'ripple_length',ripple_length);
% ncwrite(init_file,'dmix_offset',dmix_offset);
% ncwrite(init_file,'dmix_slope',dmix_slope);
% ncwrite(init_file,'dmix_time',dmix_time);
% 
% ncwrite(init_file,'plant_height',plant_height);
% ncwrite(init_file,'plant_diameter',plant_diameter);
% ncwrite(init_file,'plant_density',plant_density);
% ncwrite(init_file,'plant_thickness',plant_thickness);
% ncwrite(init_file,'marsh_mask',marsh_mask);

%close file
disp(['created ', init_file])

end
