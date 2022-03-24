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


    A_b = readField(filename, 'b_xy', k);
    A_u = readField(filename, 'u_xy', k);
    A_v = readField(filename, 'v_xy', k);
    A_w = readField(filename, 'w_xy', k);


    %A_zvort = ddx(A_v,dx) - ddy(A_u,dy);

    %A_Ri_g = ddz(A_b,dzf,1) ./  ( ddz(A_u,dzf,1).^2 + ddz(A_w,dzf,0).^2);

    A_omega_x = readField(filename, 'omegaX_xy', k);
    A_omega_y = readField(filename, 'omegaY_xy', k);
    A_omega_z = readField(filename, 'omegaZ_xy', k);

    if IC_type == 5 || IC_type == 6 || IC_type == 7 % Finite Front
      % Make full field for b
      for i = 1:size(A_b,1)
        A_b(i,:) = A_b(i,:) + dTHdx*(x(i)-Lx/2); % To get full buoyancy field
      end
    end



    %% Plot more pure fields

    dat = A_b*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,y,dat');
    colourbar('$Ro \;b$',rdbu,'Centred');
    ylabel('$y$')

    dat = A_w;
    axes(ax(2));
    pcolor_coarse(gca,x,y,dat');
    colourbar('$w$',rdbu,'Centred');

    dat = A_u; % Buoyancy Flux
    axes(ax(3));
    pcolor_coarse(gca,x,y,dat');
    colourbar('$u$',rdbu,'Centred');
    xlabel('$x$')
    ylabel('$y$')

    dat = A_v; % Vorticity
    axes(ax(4));
    pcolor_coarse(gca,x,y,dat');
    colourbar('$v$',rdbu,'Centred');
    xlabel('$x$')


    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/XY_Prim' sprintf('%04i', k) 'lres.png'],[2,2])


    %% Make figure with omega_z, PV, etc...

    dat = A_b*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,y,dat');
    colourbar('$Ro \;b$',rdbu,'Centred');
    ylabel('$y$')

    dat = A_omega_z;
    axes(ax(2));
    pcolor_coarse(gca,x,y,dat');
    colourbar('$\omega_z$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);

    dat = A_omega_x;
    axes(ax(3));
    pcolor_coarse(gca,x,y,dat');
    colourbar('$\omega_x$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')
    ylabel('$y$')

    dat = A_omega_y;
    axes(ax(4));
    pcolor_coarse(gca,x,y,dat');
    colourbar('$\omega_y$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$x$')

    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/XY_Grad' sprintf('%04i', piter) 'lres.png'],[2,2])


    close(f);


end


dT = h5readatt(filename,['/u_xy/' sprintf('%04i',Nt)],'Time') / (Ro*delta)  / double(Nt);
compile_movie('movie/XY_Prim','XY_Prim', dT);
compile_movie('movie/XY_Grad','XY_Grad', dT);
