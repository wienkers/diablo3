% Loads Movie Variables and Such

mkdir(full_dir, '/movie');

% Set the filename
filename = [full_dir '/movie.h5'];
Nt = h5readatt(filename,'/u_xy','SAMPLES') - 1;
fprintf('%i Samples Available (minus 2) \n',Nt)

if is3D
    filename_mean = [full_dir '/mean_xz.h5'];
end


%% Modify Space-time Resolution
% Figure out dt
if exist('dT_movie','var') % Set to be ~60 frames per period
    dT_movie = Ro*delta * dT_movie;
else
    dT_movie = 0.1 * (Ro*delta);
end

time1 = h5readatt(filename,['/u_xz/' sprintf('%04i',1)],'Time');
time2 = h5readatt(filename,['/u_xz/' sprintf('%04i',2)],'Time');
dt = time2 - time1;
dk = max(1,round(dT_movie/dt));

jStop = fix(Nt/dk);

t = zeros([1,jStop]);


%% Load useful colourmaps
rdbu = flipud(cbrewer('div', 'RdBu', 128));
piyg = flipud(cbrewer('div', 'PiYG', 128));
rd = cbrewer('seq', 'Reds', 128); % For Positives
    gn = cbrewer('seq','BuGn',128);
bu = flipud(cbrewer('seq', 'Blues', 128)); % For Negatives
