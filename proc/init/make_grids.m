% Makes the Grids

if exist([full_dir '/mean.mat'],'file')
    load([full_dir '/mean.mat'], 'gzf'); gzf = gzf(:);
else
    gzf = readField([full_dir '/mean.h5'], 'gzf', 1); gzf = gzf(:);   % Nz  , 1 and Nz on boundary
end

gz  = h5read([full_dir '/grid.h5'], '/grids/y');  gz  = gz(:);    % Nz+1, only 2:Nz in domain

dzf = zeros([Nz,1]);
dz = zeros([Nz,1]);


for j = 2:Nz
  dz(j) = gzf(j) - gzf(j-1); % Centred on the gz point
end

for j = 1:Nz
  dzf(j) = (gz(j+1) - gz(j)); % Centred on the gzf point
end


x = linspace(0,Lx,Nx+1); x = x(:);
x(end) = [];

dx = Lx/(Nx+1); % Last point is periodic


y = linspace(0,Ly,Ny+1); y = y(:);
y(end) = [];

dy = Ly/(Ny+1); % Last point is periodic


gz = gz(1:end-1); % To be same size as GZ data...


