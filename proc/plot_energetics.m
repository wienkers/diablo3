% Plots the Domain-Averaged Energy Budget (2D/3D)
close all; clear all;

read_mean;

[~,~,~] = mkdir(full_dir, '/budget');

colours;


%% Budget (with ddt, etc showing net budget closure)

f = figure();
% TKE
subplot(3,1,1);
plot(t_f, ddt(tke_mean,time), 'k-'); hold on;
plot(t_f, prod_mean, 'Color', CS{1});
plot(t_f, thw_mean, 'Color', CS{2});
plot(t_f, -epsilon_mean, 'Color', CS{3});
    error_TKE = ddt(tke_mean,time) - prod_mean - thw_mean + epsilon_mean;
plot(t_f, error_TKE, 'r--');
ylabel('$\mathrm{d}_t \left<E_K\right> \; [H^2 M^6/f^3]$')
xlabel('$t_f \; [f^{-1}]$')
xlim([-inf,max(t_f)]);
legend('$E_K$','$\mathcal{P}$','$\mathcal{B}$','$-\varepsilon_t$','Error')

% MKE
subplot(3,1,2);
plot(t_f, ddt(mke_mean,time), 'k-'); hold on;
plot(t_f, -prod_mean, 'Color', CS{1});
plot(t_f, buoy_m_adj_mean, 'Color', CS{2});
plot(t_f, -epsilon_m_mean, 'Color', CS{3});
    error_MKE = ddt(mke_mean,time) + prod_mean - buoy_m_adj_mean + epsilon_m_mean;
plot(t_f, error_MKE, 'r--');
ylabel('$\mathrm{d}_t \left<\overline{E_K}\right> \; [H^2 M^6/f^3]$')
xlabel('$t_f \; [f^{-1}]$')
xlim([-inf,max(t_f)]);
legend('$\overline{E_K}$','$-\mathcal{P}$','$\mathcal{B}_m^\dagger$','$-\varepsilon_m$', 'Error')

% MPE
subplot(3,1,3);
plot(t_f, ddt(pe_mean,time), 'k-'); hold on;
plot(t_f, -buoy_m_adj_mean, 'Color', CS{1});
plot(t_f, -thw_mean,'Color', CS{2});
plot(t_f, -epsilon_p_mean, 'Color', CS{3});
    error_MPE = ddt(pe_mean,time) + buoy_m_adj_mean + thw_mean + epsilon_p_mean;
plot(t_f, error_MPE, 'r--');
ylabel('$\mathrm{d}_t \left<E_P\right> \; [H^2 M^6/f^3]$')
xlabel('$t_f \; [f^{-1}]$')
xlim([-inf,max(t_f)]);
legend('$E_P$','$-\mathcal{B}_m^\dagger$','$-\mathcal{B}$','$-\varepsilon_p$', 'Error')

save_fig(f,[full_dir '/budget/KE_Budget_Closed.eps'],[2,2])




%% Energy (i.e. not showing closure)
f = figure();
% TKE
subplot(3,1,1);
plot(t_f, tke_mean, 'k-'); hold on;
plot(t_f, cumtrapz(time,prod_mean), 'Color', CS{1});
plot(t_f, cumtrapz(time,thw_mean), 'Color', CS{2});
plot(t_f, cumtrapz(time,-epsilon_mean), 'Color', CS{3});
plot(t_f, tke_mean + cumtrapz(time, - prod_mean - thw_mean + epsilon_mean), 'r--');
ylabel('$\left<E_K\right> \; [H^2 M^4/f^2]$')
xlabel('$t_f \; [f^{-1}]$')
xlim([-inf,max(t_f)]);
legend('$E_K$','$\mathcal{P}$','$\mathcal{B}$','$-\varepsilon_t$','Error')

% MKE
subplot(3,1,2);
plot(t_f, mke_mean, 'k-'); hold on;
plot(t_f, cumtrapz(time,-prod_mean), 'Color', CS{1});
plot(t_f, cumtrapz(time,buoy_m_adj_mean), 'Color', CS{2});
plot(t_f, cumtrapz(time,-epsilon_m_mean), 'Color', CS{3});
plot(t_f, mke_mean + cumtrapz(time,prod_mean - buoy_m_adj_mean + epsilon_m_mean), 'r--');
ylabel('$\left<\overline{E_K}\right> \; [H^2 M^4/f^2]$')
xlabel('$t_f \; [f^{-1}]$')
xlim([-inf,max(t_f)]);
legend('$\overline{E_K}$','$-\mathcal{P}$','$\mathcal{B}_m^\dagger$','$-\varepsilon_m$', 'Error')

% MPE
subplot(3,1,3);
plot(t_f, pe_mean, 'k-'); hold on;
plot(t_f, cumtrapz(time,-buoy_m_adj_mean), 'Color', CS{1});
plot(t_f, cumtrapz(time,-thw_mean), 'Color', CS{2});
plot(t_f, cumtrapz(time,-epsilon_p_mean), 'Color', CS{3});
plot(t_f, pe_mean + cumtrapz(time,buoy_m_adj_mean + thw_mean + epsilon_p_mean), 'r--');
ylabel('$\left<E_P\right> \; [H^2 M^4/f^2]$')
xlabel('$t_f \; [f^{-1}]$')
xlim([-inf,max(t_f)]);
legend('$E_P$','$-\mathcal{B}_m^\dagger$','$-\mathcal{B}_m$','$-\varepsilon_p$', 'Error')

save_fig(f,[full_dir '/budget/KE_Budget_Int.eps'],[2,2])



%% Production Terms

f = figure();
% TKE
plot(t_f, prod_mean, 'k-'); hold on;
plot(t_f, uu_dudx_mean, 'Color', CS{1});
plot(t_f, ww_dwdz_mean, 'Color', CS{2});
plot(t_f, vu_dvdx_mean, 'Color', CS{3});
plot(t_f, wu_dwdx_mean, 'Color', CS{4});
plot(t_f, uw_dudz_mean, 'Color', CS{5});
plot(t_f, vw_dvdz_mean, 'Color', CS{6});

ylabel('$\left<\mathcal{P}\right> \; [H^2 M^6/f^3]$')
xlabel('$t_f \; [f^{-1}]$')
xlim([-inf,max(t_f)]);
legend('$\mathcal{P}$', ...
        '$\left< \overline{u''u''}\partial_x\bar{u} \right>$',...
        '$\left< \overline{w''w''}\partial_z\bar{w} \right>$',...
        '$\left< \overline{v''u''}\partial_x\bar{v} \right>$',...
        '$\left< \overline{w''u''}\partial_x\bar{w} \right>$',...
        '$\left< \overline{u''w''}\partial_z\bar{u} \right>$',...
        '$\left< \overline{v''w''}\partial_z\bar{v} \right>$' )

save_fig(f,[full_dir '/budget/KE_Production.eps'])




%% BPE/APE

f = figure();
plot(t_f, pe_mean, 'k-'); hold on;
plot(t_f, bpe_mean, 'Color', CS{1});
plot(t_f, ape_mean, 'Color', CS{2});

ylabel('$\left<E_P\right> \; [H^2 M^6/f^3]$')
xlabel('$t_f \; [f^{-1}]$')
xlim([-inf,max(t_f)]);
legend('$E_P$', ...
        '$E_{P,b}$',...
        '$E_{P,a}$')

save_fig(f,[full_dir '/budget/PE_APE.eps'])


% f = figure();
% plot(t_f, ddt(pe_mean,time), 'k-'); hold on;
% plot(t_f, ddt(bpe_mean,time), 'Color', CS{1});
% plot(t_f, ddt(ape_mean,time), 'Color', CS{2});
% 
% ylabel('$\mathrm{d}_t \left<E_P\right> \; [H^2 M^6/f^2]$')
% xlabel('$t_f \; [f^{-1}]$')
% xlim([-inf,max(t_f)]);
% legend('$E_P$', ...
%         '$E_{P,b}$',...
%         '$E_{P,a}$')
% 
% save_fig(f,[full_dir '/budget/PE_APE_dt.eps'],[2,2])













