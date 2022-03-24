%% Plots mean TimeSeries of various quantities
%   NOTE: See old/archived file for 3D and Lagrange and x=0 line implementation
close all; clean;

xl = 80;

read_mean;

mkdir(full_dir, '/timeseries');

GZF = repmat(gzf.',[1,length(time)]);

vme_tot = vme + dVdz*(GZF-Lz/2);


%% Plot the mean TimeSeries

f = plotTimeSeries(time/(Ro*delta),gzf,ume/(pi*delta),'$\bar{u}/(\pi\delta_0)\; [\mathcal{L}_d f]$');  xlim([0,xl]);
save_fig(f,[full_dir '/timeseries/uMeanTS.png'])

f = plotTimeSeries(time/(Ro*delta),gzf,vme_tot/(pi*delta),'$\bar{v}/(\pi\delta_0)\; [\mathcal{L}_d f]$');  xlim([0,xl]);
save_fig(f,[full_dir '/timeseries/vMeanTS.png'])


vme_diff = vme_tot - repmat(vme_tot(:,1),[1,length(time)]);

Mag = sqrt(ume.^2 + vme_diff.^2);

D1 = - ume .* ddz(uw,dzf,0) ./ Mag         * Lx/(pi*delta);
D2 = - vme_diff .* ddz(wv,dzf,0) ./ Mag    * Lx/(pi*delta);


f = plotTimeSeries(time/(Ro*delta),gzf,D1 * Gamma  ,'$\mathcal{D}_u/|\bar{\mathbf{u}}|\; [\mathcal{L}_d f^2]$',true);  xlim([0,xl]);
save_fig(f,[full_dir '/timeseries/udzuwTS.png'])

f = plotTimeSeries(time/(Ro*delta),gzf,D2 * Gamma  ,'$\mathcal{D}_v/|\bar{\mathbf{u}}|\; [\mathcal{L}_d f^2]$',true);  xlim([0,xl]);
save_fig(f,[full_dir '/timeseries/vdzvwTS.png'])

f = plotTimeSeries(time/(Ro*delta),gzf,(D1+D2) * Gamma  ,'$(\mathcal{D}_u + \mathcal{D}_v)/|\bar{\mathbf{u}}|\; [\mathcal{L}_d f^2]$',true);  xlim([0,xl]);
save_fig(f,[full_dir '/timeseries/udzuvwSumTS.png'])



D1_mean = trapz(gzf,D1,1);
D2_mean = trapz(gzf,D2,1);

Mag_mean = trapz(gzf,Mag,1);

save([full_dir '/MomentumDamping.mat'],'D1_mean','D2_mean','Mag_mean','time','Ro','delta');








