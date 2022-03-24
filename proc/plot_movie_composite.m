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

    A_b_yz = readField(filename, 'b_yz', k);
    A_u_yz = readField(filename, 'u_yz', k);
    A_v_yz = readField(filename, 'v_yz', k);
    A_w_yz = readField(filename, 'w_yz', k);

    A_b_xz = readField(filename, 'b_xz', k);
    A_u_xz = readField(filename, 'u_xz', k);
    A_v_xz = readField(filename, 'v_xz', k);
    A_w_xz = readField(filename, 'w_xz', k);





    for j = 1:size(A_v_xz,2)
        A_v_xz(:,j) = A_v_xz(:,j)+dVdz*(gzf(j)-Lz/2);  % To get full velocity field
    end
    for j = 1:size(A_v_yz,2)
        A_v_yz(:,j) = A_v_yz(:,j)+dVdz*(gzf(j)-Lz/2);  % To get full velocity field
    end


    if IC_type == 5 || IC_type == 7 % Finite Front
      % Make full field for b
      for i = 1:size(A_b_xz,1)
        A_b_xz(i,:) = A_b_xz(i,:) + dTHdx*(x(i)-Lx/2); % To get full buoyancy field
      end
    elseif IC_type == 4 || IC_type == 6 % Infinite Front
      % Make perturbation field for b
      for j = 1:size(A_b_xz,2)
        A_b_xz(:,j) = A_b_xz(:,j) - Ri*(gzf(j)-Lz/2); % To get perturbation field
      end
    end


    A_omega_x_xz = readField(filename, 'omegaX_xz', k);
    A_omega_x_yz = readField(filename, 'omegaX_yz', k);
    A_omega_y_xz = readField(filename, 'omegaY_xz', k);
    A_omega_y_yz = readField(filename, 'omegaY_yz', k);



    %% Primaries

    dat = A_u_xz;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$u$',rdbu,'Centred');
    hold on; plot([Lx/2,Lx/2],[0,1],'k--'); hold off;
    ylabel('$z$')

    dat = A_u_yz;
    axes(ax(2));
    pcolor_coarse(gca,y,gzf,dat');
    colourbar('$u$',rdbu,'Centred');

    dat = A_v_xz;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$v$',rdbu,'Centred');
    hold on; plot([Lx/2,Lx/2],[0,1],'k--'); hold off;
    xlabel('$x$')
    ylabel('$z$')

    dat = A_v_yz;
    axes(ax(4));
    pcolor_coarse(gca,y,gzf,dat');
    colourbar('$v$',rdbu,'Centred');
    xlabel('$y$')

    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/XZandYZ_Prim' sprintf('%04i', piter) 'lres.png'],[2,2])



    %% Gradients

    dat = A_b_xz*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ro \;b$',rdbu,'Centred');
    hold on; plot([Lx/2,Lx/2],[0,1],'k--'); hold off;
    ylabel('$z$')

    dat = A_b_yz*Ro;
    axes(ax(2));
    pcolor_coarse(gca,y,gzf,dat');
    colourbar('$Ro \;b$',rdbu,[-1,1]);

    dat = A_omega_x_xz;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\omega_x$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    hold on; plot([Lx/2,Lx/2],[0,1],'k--'); hold off;
    xlabel('$x$')
    ylabel('$z$')

    dat = A_omega_x_yz;
    axes(ax(4));
    pcolor_coarse(gca,y,gzf,dat');
    colourbar('$\omega_x$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$y$')

    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/XZandYZ_GradX' sprintf('%04i', piter) 'lres.png'],[2,2])



    dat = A_b_xz*Ro;
    axes(ax(1));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$Ro \;b$',rdbu,'Centred');
    hold on; plot([Lx/2,Lx/2],[0,1],'k--'); hold off;
    ylabel('$z$')

    dat = A_b_yz*Ro;
    axes(ax(2));
    pcolor_coarse(gca,y,gzf,dat');
    colourbar('$Ro \;b$',rdbu,[-1,1]);

    dat = A_omega_y_xz;
    axes(ax(3));
    pcolor_coarse(gca,x,gzf,dat');
    colourbar('$\omega_y$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    hold on; plot([Lx/2,Lx/2],[0,1],'k--'); hold off;
    xlabel('$x$')
    ylabel('$z$')

    dat = A_omega_y_yz;
    axes(ax(4));
    pcolor_coarse(gca,y,gzf,dat');
    colourbar('$\omega_y$',rdbu,[-2*rms(dat(:)),2*rms(dat(:))]);
    xlabel('$y$')

    set(ax(1:2),'XTickLabel',[]);
    set(ax([2,4]),'YTickLabel',[]);
    set(ax(1:4),'layer','top');
    save_fig(f,[full_dir '/movie/XZandYZ_GradY' sprintf('%04i', piter) 'lres.png'],[2,2])


    close(f);


end


dT = h5readatt(filename,['/u_xy/' sprintf('%04i',Nt)],'Time') / (Ro*delta)  / double(Nt);
compile_movie('movie/XZandYZ_Prim','XZandYZ_Prim', dT);
compile_movie('movie/XZandYZ_GradX','XZandYZ_GradX', dT);
compile_movie('movie/XZandYZ_GradY','XZandYZ_GradY', dT);
