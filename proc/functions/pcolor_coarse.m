function ax = pcolor_coarse(ax,x,z,dat)

Nx_lim = 512; skip_x = 1;
Nz_lim = 512; skip_z = 1;

Nx = length(x);
Nz = length(z);

% Ensure x is _second_ dimension:
if size(dat,1) == Nx && Nx ~= Nz
    dat = dat.';
end

if fix(Nx/Nx_lim) >= 2
    skip_x = floor(Nx/Nx_lim);
end

if fix(Nz/Nz_lim) >= 2
    skip_z = floor(Nz/Nz_lim);
end

dat = conv2(dat, ones([skip_z, skip_x])/(skip_z*skip_x), 'same');

xg = squeeze(x(1:skip_x:end));
x_end = false;
if x(end) ~= xg(end)
    xg(end+1) = x(end);
    x_end = true;
end

zg = squeeze(z(1:skip_z:end));
z_end = false;
if z(end) ~= zg(end)
    zg(end+1) = z(end);
    z_end = true;
end

if x_end && z_end
    dat = dat([1:skip_z:end end], [1:skip_x:end end]);
elseif x_end
    dat = dat(1:skip_z:end, [1:skip_x:end end]);
elseif z_end
    dat = dat([1:skip_z:end end], 1:skip_x:end);
else
    dat = dat(1:skip_z:end, 1:skip_x:end);
end


ax = pcolor(ax, xg,zg,dat); shading interp;


end
