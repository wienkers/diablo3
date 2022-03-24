function sim = compute_critical(sim,saveout,base_dir)
%% Computes averages/integrals through tau_crit

if nargin < 3
    clear base_dir;
end

plotDiagnostics = true;

read_mean;


%% Setup the storage

try
    sim.Ro;
catch

    sim.Ro = [];

    sim.is3D = [];

    sim.sigma = [];

    sim.tau_crit = [];

    sim.Upar_crit = [];
    sim.Upar_mean = [];

    sim.TKE_crit = [];
    sim.TKE_mean = [];

    sim.intProd_crit = [];
    sim.intBuoy_crit = [];
    sim.intProd_10 = [];
    sim.intBuoy_10 = [];

    sim.ProdOver2E = [];
    sim.BuoyOver2E = [];

    sim.intdz2wb_crit = [];
    sim.intdzUoRo_crit = [];
    sim.intdz2wb_10 = [];
    sim.intdzUoRo_10 = [];

    sim.intdzvw_crit = [];
    sim.intdzuw_crit = [];

    sim.intdV_min = [];

end

%load('./critical.mat')


sim.Ro = [sim.Ro Ro];

sim.is3D = [sim.is3D is3D];


%% Compute k_threshold above which SI is all linearly stable
%fun = @(k) abs(real(computeBoundedDispersionSemiViscousNHS(k,Re,Ro,1,1,0,true)));

%k_high = fminbnd(fun,10,1000)


%% Compute tau_crit

%k_high = 300;
%E_thresh = 1e-7;
slope0_thresh = 0;
diss_crit = 2e-4;


% kx = (0:Nx/2-1) * 2*pi/Lx;
% KX = repmat(kx.',[1,Nt]);
%
% FTx_uu_high = FTx_uu;
% FTx_uu_high(KX < k_high) = 0;
%
% highKxEnergy = trapz(kx,FTx_uu_high,1);
% totalEnergy = trapz(kx,FTx_uu,1);

% Start time, at the bottom of the dip
tau_ = intersections(time,ddt(epsilon_int,time),time,ones(size(time))*slope0_thresh,false);
if tau_(1) < 1e-10
    tau_(1) = [];
end
tau_start = tau_(1)
[~,ind_start] = min(abs(tau_start - time));


% tau_crit = intersections(time,highKxEnergy,time,ones(size(time))*E_thresh,false);
% if tau_crit(1) < 5
%     tau_crit(1) = [];
% end
%tau_crit = tau_(2)
%[~,ind_crit] = min(abs(tau_crit - time));

tau_crit = intersections(time,epsilon_int,time,ones(size(time))*diss_crit,false);
if tau_crit(1) < 1
    tau_crit(1) = [];
end
tau_crit = tau_crit(1)
sim.tau_crit = [tau_crit sim.tau_crit];
[~,ind_crit] = min(abs(tau_crit - time));


% 10 periods after tau_crit
tau_10 = (tau_crit + 2*pi*Ro*delta*10);
[~,ind_10] = min(abs(tau_10 - time));



% Time when mid-way up the linear growth
tau_mid = tau_start + 0.5*(tau_crit - tau_start);
[~,ind_mid] = min(abs(tau_mid - time));

mask_linear = time > tau_mid & time < tau_crit;

if plotDiagnostics
    figure()
    semilogy(time,tke_int); hold on
    %semilogy(time,totalEnergy)
    %plot(time,highKxEnergy./totalEnergy)
    %semilogy(time,ddt(highKxEnergy,time))
    plot(tau_crit,tke_int(ind_crit),'x')
    plot(tau_mid,tke_int(ind_mid),'x')
    plot(tau_start,tke_int(ind_start),'x')
    ylabel('Energy')
end


%% Compute sigma averaged on (0, tau_crit)

sigma = (log(tke_int(ind_crit)) - log(tke_int(ind_mid))) / (tau_crit - tau_mid);
sigma = sigma/2

sim.sigma = [sim.sigma sigma];

if plotDiagnostics
    plot(time(mask_linear),tke_int(ind_crit) * exp(2*sigma*(time(mask_linear) - time(ind_crit))),'r--')
end



%% Compute Normalised GS Shear and Buoyancy Production (in linear regime)

Po2E = -trapz(gzf,wv,1)./(2*tke_int);
Bo2E = thw_int./(2*tke_int);


Po2E_mean = mean(Po2E(mask_linear))
sim.ProdOver2E = [sim.ProdOver2E Po2E_mean];
Bo2E_mean = mean(Bo2E(mask_linear))
sim.BuoyOver2E = [sim.BuoyOver2E Bo2E_mean];



if plotDiagnostics
    figure()
    plot(time,Po2E); hold on;
    plot(time,Bo2E);
    plot(tau_crit,Po2E(ind_crit),'x')
    plot(tau_mid,Po2E(ind_mid),'x')
    ylabel('Po2E and Bo2E for mean')
    legend('Po2E','Bo2E')
end




%% Compute U_||,rms & TKE

% After 10 periods: tau_crit < time < tau_crit + 2*pi*Ro*10

mask_nonlinear = time > tau_crit & time < tau_10;


U_parRMS = sqrt(urms_int.^2 + wrms_int.^2);

if plotDiagnostics
    figure()
    semilogy(time,U_parRMS); hold on;
    plot(tau_crit,U_parRMS(ind_crit),'x')
    plot(tau_mid,U_parRMS(ind_mid),'x')
    ylabel('U_par')
end


U_parRMS_atCrit = U_parRMS(ind_crit);

U_par_atCrit = sqrt(2) * U_parRMS_atCrit
sim.Upar_crit = [sim.Upar_crit U_par_atCrit];


U_par_mean = sqrt(2) * mean(U_parRMS(mask_nonlinear))
sim.Upar_mean = [sim.Upar_mean U_par_mean];


TKE_atCrit = tke_int(ind_crit)
sim.TKE_crit = [sim.TKE_crit TKE_atCrit];

TKE_mean = mean(tke_int(mask_nonlinear))
sim.TKE_mean = [sim.TKE_mean TKE_mean];


%% Compute Integrated KE Budget Terms

Prod = -trapz(gzf,wv,1);
Buoy = thw_int;

if plotDiagnostics
    figure()
    semilogy(time,Prod); hold on;
    plot(time,Buoy)
    plot(tau_crit,Prod(ind_crit),'x')
    legend('Prod','Buoy')
end

P_toCrit = trapz(time(1:ind_crit),Prod(1:ind_crit))
sim.intProd_crit = [sim.intProd_crit P_toCrit];

B_toCrit = trapz(time(1:ind_crit),Buoy(1:ind_crit))
sim.intBuoy_crit = [sim.intBuoy_crit B_toCrit];

P_int = trapz(time(1:ind_10),Prod(1:ind_10))
sim.intProd_10 = [sim.intProd_10 P_int];

B_int = trapz(time(1:ind_10),Buoy(1:ind_10))
sim.intBuoy_10 = [sim.intBuoy_10 B_int];


%% Compute the Buoyancy Budget Term dz2wb

dz2wb = trapz(gzf,abs(d2dz(thw,dzf,0)));

dzUoRo = trapz(gzf,abs(ddz(ume,dzf,1)))/(Ro*delta);

if plotDiagnostics
    figure()
    semilogy(time,dz2wb); hold on
    plot(time,dzUoRo)
    legend('dz2wb','dzU0Ro')
end

dz2wb_toCrit = trapz(time(1:ind_crit),dz2wb(1:ind_crit))
sim.intdz2wb_crit = [sim.intdz2wb_crit dz2wb_toCrit];

dzUoRo_toCrit = trapz(time(1:ind_crit),dzUoRo(1:ind_crit))
sim.intdzUoRo_crit = [sim.intdzUoRo_crit dzUoRo_toCrit];

dz2wb_10 = trapz(time(1:ind_10),dz2wb(1:ind_10))
sim.intdz2wb_10 = [sim.intdz2wb_10 dz2wb_10];

dzUoRo_10 = trapz(time(1:ind_10),dzUoRo(1:ind_10))
sim.intdzUoRo_10 = [sim.intdzUoRo_10 dzUoRo_10];


%% Compute Momentum Budget Terms

dzvw = trapz(gzf,abs(ddz(wv,dzf,0)));

dzuw = trapz(gzf,abs(ddz(uw,dzf,0)));

if plotDiagnostics
    figure()
    semilogy(time,dzvw); hold on
    plot(time,dzuw);
    legend('dzvw','dzuw')
end

dzvw_toCrit = sqrt(2) * trapz(time(1:ind_crit),dzvw(1:ind_crit))
sim.intdzvw_crit = [sim.intdzvw_crit dzvw_toCrit];

dzuw_toCrit = sqrt(2) * trapz(time(1:ind_crit),dzuw(1:ind_crit))
sim.intdzuw_crit = [sim.intdzuw_crit dzuw_toCrit];


intdV = vme(end,:) - vme(1,:);
sim.intdV_min = [sim.intdV_min min(intdV(time < tau_crit + 2*pi*Ro*delta))];

if plotDiagnostics
    figure()
    plot(time,intdV); hold on
    legend('intdV')
end




if saveout
    save([full_dir, '/critical.mat'],'sim')
end


end
