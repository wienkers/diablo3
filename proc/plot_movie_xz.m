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


    A_b = readField(filename, 'b_xz', k);
    A_u = readField(filename, 'u_xz', k);
    A_v = readField(filename, 'v_xz', k);
    A_w = readField(filename, 'w_xz', k);

    A_omega_x = readField(filename, 'omegaX_xz', k);
    A_omega_y = readField(filename, 'omegaY_xz', k);
    A_omega_z = readField(filename, 'omegaZ_xz', k);


%     A_yvort = ddz(A_u,dzf,1) - ddx(A_w,dx);
%     A_zvort = ddx(A_v,dx);
%     A_xvort = -ddz(A_v,dzf,1) - dVdz;

    %A_Ri_g = ddz(A_b,dzf,1) ./  ( ddz(A_u,dzf,1).^2 + ddz(A_w,dzf,0).^2);

    %q = PV(A_v,A_b,Ro*delta,dTHdx,dVdx,dVdz,dx,dzf,gzf);


    for j = 1:size(A_v,2)
        A_v(:,j) = A_v(:,j)+dVdz*(gzf(j)-Lz/2);  % To get full velocity field
    end

    if IC_type == 5 || IC_type == 6 || IC_type == 7 % Finite Front
      % Make full field for b
      for i = 1:size(A_b,1)
        A_b(i,:) = A_b(i,:) + dTHdx*(x(i)-Lx/2); % To get full buoyancy field
      end
    elseif IC_type == 4 % Infinite Front
      % Make perturbation field for b
      for j = 1:size(A_b,2)
        A_b(:,j) = A_b(:,j) - Ri*(gzf(j)-Lz/2); % To get perturbation field
      end
    end

    ke = A_u.^2 + A_v.^2 + A_w.^2;


    %% Plot more pure fields

    dat = A_b*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ro \;b$',rdbu,'Centred');
    ylabel('$z$')

    dat = A_w;
    axes(ax(2));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$w$',rdbu,'Centred');

    dat = A_u;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$u$',rdbu,'Centred');
    xlabel('$x$')
    ylabel('$z$')

    dat = A_v;
    axes(ax(4));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$v$',rdbu,'Centred');
    xlabel('$x$')


    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/XZ_Prim' sprintf('%04i', piter) 'lres.png'],[2,2])



    %% Make figure with omega_z, PV, etc...

    dat = A_b*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ro \;b$',rdbu,'Centred');
    ylabel('$z$')

    dat = A_omega_z;
    axes(ax(2));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\omega_z$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);

    dat = A_omega_x;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\omega_x$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')
    ylabel('$z$')

    dat = A_omega_y;
    axes(ax(4));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\omega_y$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')

    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/XZ_Grad' sprintf('%04i', piter) 'lres.png'],[2,2])


    close(f);


end


dT = h5readatt(filename,['/u_xy/' sprintf('%04i',Nt)],'Time') / (Ro*delta)  / double(Nt);
compile_movie('movie/XZ_Grad','XZ_Grad', dT);
compile_movie('movie/XZ_Prim','XZ_Prim', dT);
