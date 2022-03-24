close all; clear all;

set_inputs;
make_grids;
dT_movie = 0.05;
init_movie;


tstart = 0  *  (Ro*delta);
parfor piter = 1:jStop
    k = 1 + (piter-1)*dk;

    t(piter) = read_time(filename,k);

    f = figure();
    ax = subplots(2,2,[0.015 0.035],[0.085,0.1],[0.06,0.05]);
    figtitle = sgtitle(''); figtitle.FontSize = 12;
    figtitle.String = sprintf('$t = %.2f \\; [1/f]$',t(piter)/(Ro*delta));

    if t(piter) < tstart - 1e-4
        close(f);
        continue;
    end


    %M_q = readField(filename_mean, 'q_xz', k); % NO LONGER SAVED
    M_b = readField(filename_mean, 'bme_xz', k);

    M_u = readField(filename_mean, 'ume_xz', k);
    M_v = readField(filename_mean, 'vme_xz', k);
    M_w = readField(filename_mean, 'wme_xz', k);

    M_uu = readField(filename_mean, 'uu_xz', k);
    M_vv = readField(filename_mean, 'vv_xz', k);
    M_ww = readField(filename_mean, 'ww_xz', k);
    
    M_chi = readField(filename_mean, 'chi_xz', k);


    [q_m, q_BC_m, q_vort_m] = PV(M_v,M_b,Gamma,dTHdx,dVdx,dVdz,dx,dzf,gzf);


    for j = 1:size(M_v,2)
        M_v(:,j) = M_v(:,j)+dVdz*(gzf(j)-Lz/2);  % To get full velocity field
    end

    if IC_type == 5 || IC_type == 6 || IC_type == 7 % Finite Front
      % Make full field for b
      for i = 1:size(M_b,1)
        M_b(i,:) = M_b(i,:) + dTHdx*(x(i)-Lx/2); % To get full buoyancy field
      end
    elseif IC_type == 4 % Infinite Front
      % Make perturbation field for b
      for j = 1:size(M_b,2)
        M_b(:,j) = M_b(:,j) - Ri*(gzf(j)-Lz/2); % To get perturbation field
      end
    end


    %% Secondary Quantities

    shear_m = sqrt( ddz(M_v, dzf, 'Neumann').^2 + ddz(M_u, dzf, 'Neumann').^2);
    Ri_m = ddz(M_b, dzf, 'Neumann') ./ shear_m;



    %% TKE Quantities

    tke = 0.5 * (M_uu + M_vv + M_ww);

    M_epsilon = readField(filename_mean, 'epsilon_xz', k);

    M_buoy = readField(filename_mean, 'thw_xz', k);

    % Compute Production
    uw = readField(filename_mean, 'uw_xz', k);
    wv = readField(filename_mean, 'wv_xz', k);
    uv = readField(filename_mean, 'uv_xz', k);

    stp = -M_uu .* ddx(M_u,dx)         +  -M_ww .* ddz(M_w, dzf, 'Dirichlet');
    lsp = -uv .* ddx(M_v,dx)                +  -uw .* ddx(M_w,dx);
    vsp = -uw .* ddz(M_u, dzf, 'Neumann')   +  -wv .* ddz(M_v, dzf, 'Neumann');

    prod = stp + lsp + vsp;


    %% MKE Quantities

    mke = 0.5 * (M_u.^2 + M_v.^2 + M_w.^2);

    try
        M_epsilon_m = readField(filename_mean, 'epsilon_m_xz', k);
    catch
        M_epsilon_m = zeros(size(tke));
    end

    buoy_m = M_b .* M_w;


    %% Plot Basic Overview

    dat = M_b*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ro \;\bar{b}$',rdbu,'Centred');
    ylabel('$z$')

    dat = Ri_m;
    axes(ax(2));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ri_m$',rdbu,[-1.5, 1.5]);

    dat = M_u;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\bar{u}$',rdbu,'Centred');
    xlabel('$x$')
    ylabel('$z$')

    dat = M_v;
    axes(ax(4));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\bar{v}$',rdbu, 'Centred');
    xlabel('$x$')


    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/meanXZ_Summary' sprintf('%04i', k) 'lres.png'],[2,2])



    %% Plot PV Summary

    dat = M_b*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ro \;\bar{b}$',rdbu,'Centred');
    ylabel('$z$')

    dat = q_BC_m;
    axes(ax(2));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\tilde{q}_b$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);

    dat = q_vort_m;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\tilde{q}_\omega$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')
    ylabel('$z$')

    dat = q_m; % PV of the _mean_ !
    axes(ax(4));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\tilde{q}$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')


    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/meanXZ_PV' sprintf('%04i', k) 'lres.png'],[2,2])



    %% Plot TKE Summary

    dat = tke;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$E_K$',rd,[0,2*rms(dat(:))]);
    ylabel('$z$')

    dat = prod;
    axes(ax(2));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\mathcal{P}$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);

    dat = M_buoy;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\mathcal{B}$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')
    ylabel('$z$')

    dat = M_epsilon;
    axes(ax(4));
    pcolor_coarse(gca,x,gzf,dat');
    %set(ax(4),'colorscale','log');
    colourbar('$\varepsilon_t$',rdbu,'Percentile98')% rd,[1e-4*rms(dat(:)),10*rms(dat(:))]);
    xlabel('$x$')


    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/meanXZ_TKE' sprintf('%04i', k) 'lres.png'],[2,2])


    %% Plot TKE Production

    dat = tke;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$E_K$',rd,[0,2*rms(dat(:))]);
    ylabel('$z$')

    dat = stp;
    axes(ax(2));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\mathcal{P}_{st}$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);

    dat = lsp;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\mathcal{P}_{ls}$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')
    ylabel('$z$')

    dat = vsp;
    axes(ax(4));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\mathcal{P}_{vs}$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')


    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/meanXZ_TKE_Prod' sprintf('%04i', k) 'lres.png'],[2,2])



    %% Plot MKE Summary

    dat = mke;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\overline{E}_K$',rd,'Positive');
    ylabel('$z$')

    dat = -prod;
    axes(ax(2));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$-\mathcal{P}$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);

    dat = buoy_m;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\mathcal{B}_m$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')
    ylabel('$z$')

    dat = M_epsilon_m;
    axes(ax(4));
    pcolor_coarse(gca,x,gzf,dat');
    set(ax(4),'colorscale','log');
    colourbar('$\varepsilon_m$',rd,[1e-4*rms(dat(:)),10*rms(dat(:))]);
    xlabel('$x$')


    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/meanXZ_MKE' sprintf('%04i', k) 'lres.png'],[2,2])
    
    
    
    %% Dissipation Summary

    dat = M_b*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ro \;\bar{b}$',rdbu,'Centred');
    ylabel('$z$')

    dat = Ri_m;
    axes(ax(2));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ri_m$',rdbu,[-1.5, 1.5]);

    dat = M_chi;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    set(ax(3),'colorscale','log');
    colourbar('$\chi$',rd,[1e-2*rms(dat(:)),10*rms(dat(:))]);
    xlabel('$x$')
    ylabel('$z$')

    dat = M_epsilon;
    axes(ax(4));
    pcolor_coarse(gca,x,gzf,dat');
    set(ax(4),'colorscale','log');
    colourbar('$\varepsilon_t$',rd,[1e-4*rms(dat(:)),10*rms(dat(:))]);
    xlabel('$x$')



    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/meanXZ_Diss' sprintf('%04i', k) 'lres.png'],[2,2])




    close(f);


end


dT = h5readatt(filename_mean,['/ume_xz/' sprintf('%04i',Nt)],'Time') / (Ro*delta)  / double(Nt);
compile_movie('movie/meanXZ_Summary','meanXZ_Summary', dT);
compile_movie('movie/meanXZ_PV','meanXZ_PV', dT);
compile_movie('movie/meanXZ_TKE','meanXZ_TKE', dT);
compile_movie('movie/meanXZ_TKE_Prod','meanXZ_TKE_Prod', dT);
compile_movie('movie/meanXZ_MKE','meanXZ_MKE', dT);
compile_movie('movie/meanXZ_Diss','meanXZ_Diss', dT);
