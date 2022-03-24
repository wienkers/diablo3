%% Reads in all of the Statistics in mean.h5

% Read inputs / Make Grid
set_inputs;
make_grids;


if exist([full_dir '/mean.mat'],'file')
    tau_crit_new = tau_crit;
    load([full_dir '/mean.mat']);
    % Check if tau_crit has been changed...
    if tau_crit ~= tau_crit_new
        % Update the time shift
        time = time + tau_crit;
        tau_crit = tau_crit_new;
        
        time = time - tau_crit; % By default, use offset time
        t_f = time / Gamma; % Inertial units
        t_p = t_f / (2*pi); % Inertial periods
    end
     
     return
end


filename = [full_dir '/mean.h5'];
Nt = h5readatt(filename,'/ume','SAMPLES') - 1;
fprintf('%i Samples Available (minus 1) \n', Nt)



%% Read in All Primary Quantities

time = zeros(1,Nt);
ume = zeros(Nz,Nt); vme = zeros(Nz,Nt); wme = zeros(Nz,Nt);
mke = zeros(Nz,Nt);
urms = zeros(Nz,Nt); vrms = zeros(Nz,Nt); wrms = zeros(Nz,Nt);
uv = zeros(Nz,Nt); uw = zeros(Nz,Nt); wv = zeros(Nz,Nt);
uu_dudx = zeros(Nz,Nt); ww_dwdz = zeros(Nz,Nt); vu_dvdx = zeros(Nz,Nt);
wu_dwdx = zeros(Nz,Nt); uw_dudz = zeros(Nz,Nt); vw_dvdz = zeros(Nz,Nt);
dudz = zeros(Nz,Nt); dvdz = zeros(Nz,Nt); cp = zeros(Nz,Nt); shear = zeros(Nz,Nt);
omega_x = zeros(Nz,Nt); omega_y = zeros(Nz,Nt); omega_z = zeros(Nz,Nt);
thme = zeros(Nz,Nt,N_th); dthdz = zeros(Nz,Nt,N_th); thrms = zeros(Nz,Nt,N_th);
thw = zeros(Nz,Nt,N_th); thw_m = zeros(Nz,Nt,N_th); chi_m = zeros(Nz,Nt,N_th);
th_int = zeros(1,Nt);

z_star = zeros(Nx,Nt);
dTHdz_star = zeros(Nx,Nt);
bpe_pdf = zeros(Nx,Nt);

chi_star = zeros(Nx,Nt);
shear2_star = zeros(Nx,Nt);
thw_star = zeros(Nx,Nt);
epsilon_star = zeros(Nx,Nt);
TKE_star = zeros(Nx,Nt);
KE_star = zeros(Nx,Nt);

kappa_t_rel = zeros(Nx,Nt);
phi_d = zeros(Nx,Nt);
kappa_t_rel_mean = zeros(1,Nt);
phi_d_mean = zeros(1,Nt);

bpe_mean = zeros(1,Nt); dbdz_s = zeros(1,Nt);
uz_x0 = zeros(1,Nt);  %th_min = zeros(1,Nt); th_max = zeros(1,Nt);
pe_mean = zeros(1,Nt);

epsilon = zeros(Nz,Nt); epsilon_m = zeros(Nz,Nt);
FTx_uu = zeros(Nx/2,Nt);

nu_sgs = zeros(Nz,Nt); epsilon_sgs = zeros(Nz,Nt);



th_bin = readField(filename, 'th1bin', 2);
th_bin = [th_bin; 2*th_bin(end)-th_bin(end-1)]; % These are constant in time!

for k = 1:Nt
    
    time(k) = read_time(filename,k);
    
    ume(:,k) = readField(filename, 'ume', k);           % Mean Velocities
    vme(:,k) = readField(filename, 'vme', k) + dVdz*(gz-Lz/2);   % On GZ    % Add the Background TWS!!!
    wme(:,k) = readField(filename, 'wme', k);
    
    mke(:,k) = readField(filename, 'mke', k);   % On GZ
    
    urms(:,k) = readField(filename, 'urms', k);         % RMS Velocities
    vrms(:,k) = readField(filename, 'vrms', k); % On GZ
    wrms(:,k) = readField(filename, 'wrms', k);
    
    uv(:,k) = readField(filename, 'uv', k);             % Reynolds Stresses
    uw(:,k) = readField(filename, 'uw', k);     % On GZ
    wv(:,k) = readField(filename, 'wv', k);     % On GZ
    
    uu_dudx(:,k) = readField(filename, 'uu_dudx', k);   % TKE Production
    ww_dwdz(:,k) = readField(filename, 'ww_dwdz', k); % On GZ
    vu_dvdx(:,k) = readField(filename, 'vu_dvdx', k);
    wu_dwdx(:,k) = readField(filename, 'wu_dwdx', k); % On GZ
    uw_dudz(:,k) = readField(filename, 'uw_dudz', k); % On GZ
    vw_dvdz(:,k) = readField(filename, 'vw_dvdz', k); % On GZ
    
    dudz(:,k) = readField(filename, 'dudz', k);         % Mean Velocity Gradients
    dvdz(:,k) = readField(filename, 'dvdz', k);
   
    cp(:,k) = readField(filename, 'cp', k);             % Mean Pressure
    shear(:,k) = readField(filename, 'shear', k);       % Mean Square Shear
    
    omega_x(:,k) = readField(filename, 'omega_x', k);   % RMS Vorticity
    omega_y(:,k) = readField(filename, 'omega_y', k);
    omega_z(:,k) = readField(filename, 'omega_z', k);
    
    
    
    if N_th ~= 1
        error('Warning: Not implemented with N_th > 1')
    end
    n = 1;
    thme(:,k) = readField(filename, ['thme' num2str(n,'%2.2d')], k);      % Mean Buoyancy
        thme(:,k) = thme(:,k) + mean(dTHdx*(x-Lx/2)); % Account for non-centred x-grid!
    dthdz(:,k) = readField(filename, ['dthdz' num2str(n,'%2.2d')], k);    % Mean Buoyancy Gradients
    thrms(:,k) = readField(filename, ['thrms' num2str(n,'%2.2d')], k);    % RMS Buoyancy
    thw(:,k) = readField(filename, ['thw' num2str(n,'%2.2d')], k);        % Buoyancy Production !!!
    thw_m(:,k) = readField(filename, ['thw' num2str(n,'%2.2d') '_m'], k); % Mean Buoyancy Production !!!
    % NOTE: THIS IS NOT PE DISSIPATION!
        chi_m(:,k) = readField(filename, ['pe_diss' num2str(n,'%2.2d')], k);
    
    th_int(k) = trapz(gzf, thme(:,k), 1) - mean(dTHdx*(x-Lx/2));

    pe_mean(k) = trapz(gzf,-(gzf-Lz/2).*thme(:,k)) / Lz; % NOTE: the BG mean X-gradient averages out at each Z
    bpe_pdf(:,k) = readField(filename, 'th1PDF', k);    % Quantities to reconstruct BPE PDF
        bpe_pdf(:,k) = bpe_pdf(:,k) / sum(bpe_pdf(:,k));
        % NOTE: This is the distribution function !NOT! density function
        % (i.e. not scaled by 1/dTH...)
    
    %bpe_mean(k) = readField(filename, 'BPE', k);        % Mean BPE (i.e. not summed)
    % NOTE: Run-time BPE is worse...
    [bpe_mean(k), ~, dbdz_s(k)] = compute_BPE(bpe_pdf(:,k), th_bin,  0,0, th_int(k),dTHdx,Lx,Lz);
    
    % Z Star Binned Quantities
    chi_star(:,k) = max(0,filter_zstar_field(readField(filename, 'chi_zstar', k)));    % Quantities to reconstruct turbulent diffusivity
    shear2_star(:,k) = filter_zstar_field(readField(filename, 'shear2_zstar', k));
    thw_star(:,k) = filter_zstar_field(readField(filename, 'thw_zstar', k));
    epsilon_star(:,k) = filter_zstar_field(readField(filename, 'epsilon_zstar', k));
    TKE_star(:,k) = filter_zstar_field(readField(filename, 'TKE_zstar', k));
    KE_star(:,k) = filter_zstar_field(readField(filename, 'KE_zstar', k));
    
    
    [kappa_t_rel(:,k), kappa_t_rel_mean(k), phi_d(:,k), phi_d_mean(k), z_star(:,k), dTHdz_star(:,k)] = ...
                                compute_KappaT(bpe_pdf(:,k), chi_star(:,k), th_bin, th_int(k), dTHdx,Lx,Lz, Re*Pr);
        % Mean quantities computed by finding phi_d(z*) and then averaging... (Should agree with global definition)
    
    
    
    uz_x0(k) = readField(filename, 'u1z_x0', k);        % <u1*z> at x edge (For Pressure work in BPE Equation)
    
    epsilon(:,k) = readField(filename, 'epsilon', k); % On GZ       % Dissipation Rate
    try
        epsilon_m(:,k) = readField(filename, 'epsilon_m', k); % On GZ   % Dissipation Rate (of MKE)
    catch
        epsilon_m(:,k) = zeros(size(epsilon(:,k)));
    end
    
    FTx_uu(:,k) = readField(filename, 'FTx_uu', k);     % Spatial X Spectrum


    % LES Output
    if LES
        nu_sgs(:,k) = readField(filename, 'nu_sgs', k);
        eps_sgs1 = -1 * readField(filename, 'eps_sgs1', k); % On GZ
        eps_sgs2 = -1 * readField(filename, 'eps_sgs2', k); % On GZ
            epsilon_sgs(:,k) = eps_sgs1 + eps_sgs2;
    else
        nu_sgs(:,k) = zeros(Nz,1);
        epsilon_sgs(:,k) = zeros(Nz,1);
    end


end

% Zero out first point on the GZ grid
vme(1,:) = 0;
mke(1,:) = 0;
vrms(1,:) = 0;
uw(1,:) = 0;
wv(1,:) = 0;
ww_dwdz(1,:) = 0;
wu_dwdx(1,:) = 0;
uw_dudz(1,:) = 0;
vw_dvdz(1,:) = 0;
epsilon(1,:) = 0;
epsilon_m(1,:) = 0;
epsilon_sgs(1,:) = 0;





%% Compute Secondary Quantities
% See previous versions for other quantities...

tke = 0.5*(urms.^2 + vrms.^2 + wrms.^2);

time = time - tau_crit; % By default, use offset time
t_f = time / Gamma; % Inertial units
t_p = t_f / (2*pi); % Inertial periods



%% Calculate Domain-Integrated _Mean_ (Budget) Quantities

urms_mean = sqrt(trapz(gzf,urms.^2,1) / Lz);
vrms_mean = sqrt(trapz(gzf,vrms.^2,1) / Lz);
wrms_mean = sqrt(trapz(gzf,wrms.^2,1) / Lz);
thrms_mean = sqrt(trapz(gzf,thrms.^2,1) / Lz);

hke_mean = 0.5 * trapz(gzf, (urms.^2 + vrms.^2), 1) / Lz;
vke_mean = 0.5 * dz.' * wrms.^2 / Lz; % Integrate conservatively (on GY)
tke_mean = vke_mean + hke_mean;

mke_mean = dz.' * mke / Lz; % Integrate conservatively (on GY)

uu_dudx_mean = trapz(gzf, -uu_dudx, 1) / Lz;
vu_dvdx_mean = trapz(gzf, -vu_dvdx, 1) / Lz;
ww_dwdz_mean = dz.' * -ww_dwdz / Lz; % Integrate conservatively (on GY)
wu_dwdx_mean = dz.' * -wu_dwdx / Lz; % Integrate conservatively (on GY)
uw_dudz_mean = dz.' * -uw_dudz / Lz; % Integrate conservatively (on GY)
vw_dvdz_mean = dz.' * -vw_dvdz / Lz; % Integrate conservatively (on GY)

stp_mean = uu_dudx_mean + ww_dwdz_mean;
lsp_mean = vu_dvdx_mean + wu_dwdx_mean;
vsp_mean = uw_dudz_mean + vw_dvdz_mean;

prod_mean = stp_mean + lsp_mean + vsp_mean;

thw_mean = trapz(gzf, thw, 1) / Lz;
thw_m_mean = trapz(gzf, thw_m, 1) / Lz;

buoy_m_adj_mean = thw_m_mean + dTHdx*uz_x0;  % Account for the x0 boundary term.


epsilon_mean = dz.' * epsilon / Lz; % Integrate conservatively (on GY)
epsilon_m_mean = dz.' * epsilon_m / Lz; % Integrate conservatively (on GY)

epsilon_sgs_mean = dz.' * epsilon_sgs / Lz; % Integrate conservatively (on GY)
nu_sgs_mean = trapz(gzf,nu_sgs,1);

epsilon_p_mean = 1/(Re*Pr) * ( - (thme(end,:)-thme(1,:))/Lz);
    % When dthdz = 0 at the boundaries...
    % Also, don't include dTHdx... (Maybe because that gradient isn't recognised in diablo?

ape_mean = pe_mean - bpe_mean;
phi_d_mean_bpe = ddt(bpe_mean,time); % Diapycnal Mixing Rate (Should be > 0)





%% Clean Up

clear recycle j k n kstart kend linenum f_input timename varname C att_info ans

%mean_variables = who; sprintf('''%s'',',mean_variables{:})
save([full_dir '/mean.mat'],'FTx_uu','Gamma','IC_type','KE_star','LES','Lx','Ly','Lz','N_th','Nt','Nx','Ny','Nz','Pr','Re','Re_v','Ri','Ro','TKE_star','ape_mean','base_dir','beta','bpe_mean','bpe_pdf','buoy_m_adj_mean','chi_star','cp','dTHdx','dTHdy','dTHdz_star','dVdx','dVdz','dbdz_s','delta','dthdz','dudz','dvdz','dx','dy','dz','dzf','epsilon','epsilon_m','epsilon_m_mean','epsilon_mean','epsilon_p_mean','epsilon_sgs','epsilon_sgs_mean','epsilon_star','filename','full_dir','chi_m','gz','gzf','hke_mean','is3D','kappa','kappa_t_rel','kappa_t_rel_mean','lsp_mean','mke','mke_mean','nu','nu_sgs','nu_sgs_mean','nu_v','omega_x','omega_y','omega_z','pe_mean','perturbation','phi_d','phi_d_mean','phi_d_mean_bpe','prod_mean','root_dir','shear','shear2_star','stp_mean','t_f','t_p','tau_crit','th_bin','th_int','thme','thrms','thrms_mean','thw','thw_m','thw_m_mean','thw_mean','thw_star','time','tke','tke_mean','ume','urms','urms_mean','uu_dudx','uu_dudx_mean','uv','uw','uw_dudz','uw_dudz_mean','uz_x0','vke_mean','vme','vrms','vrms_mean','vsp_mean','vu_dvdx','vu_dvdx_mean','vw_dvdz','vw_dvdz_mean','wme','wrms','wrms_mean','wu_dwdx','wu_dwdx_mean','wv','ww_dwdz','ww_dwdz_mean','x','y','z_star');




function field = filter_zstar_field(field)
%     mask = isnan(field);
%     tx    = 1:numel(field);
%     field(mask) = interp1(tx(~mask), field(~mask), tx(mask));
%     field(isnan(field)) = 1e-100;
    
    %field(field > 1e50) = 1e-100; 
    %field(field < 1e-100) = 1e-100;
end



