% Sets input parameters


% Read in the input directory
if exist('base_dir','var') % We've overridden the base_dir in script
    % Do nothing...
elseif ~exist('./post_dir','file') % base_dir
    root_dir = '';
    base_dir = pwd;
else % post_dir contains the base_dir and then if non-default, the root_dir
    fid = fopen('post_dir');
    base_dir = fgetl(fid);
    temp = fgetl(fid);
    if temp(1) ~= -1
        root_dir = temp;
    end
end

% Check if already have the correct data loaded
recycle = false;
if exist('full_dir','var') && contains(full_dir,base_dir) && exist('urms','var')
    recycle = true
    return
end

if ~exist('root_dir','var')
    root_dir = '/Users/Wienkers/Documents/Education/University_of_Cambridge/Research/Data/';
end
full_dir = [root_dir base_dir];


%% Read from input files

f_input = fopen([full_dir '/input.dat']);

linenum = 7;
C = textscan(f_input,'%s',2, 'headerlines',linenum-1);
LES = C{1}(2);
LES = strcmp(LES,'.TRUE.');

C = textscan(f_input,'%f',5, 'headerlines',2); % 2 lines below...
Re = C{1}(1);  nu = 1/Re;
beta = C{1}(2);
Lx = C{1}(3);
Ly = C{1}(4);
Lz = C{1}(5);

C = textscan(f_input,'%f',1, 'headerlines',2); % 2 lines below...
nu_v = C{1}(1) * nu;  Re_v = 1/nu_v;

C = textscan(f_input,'%f',2, 'headerlines',13);
Ri = C{1}(1);
Pr = C{1}(2);  kappa = nu/Pr;


f_input = fopen([full_dir '/input_chan.dat']);

linenum = 11;
C = textscan(f_input,'%f',2, 'headerlines',linenum-1);
IC_type = C{1}(1);
perturbation = C{1}(2);

C = textscan(f_input,'%f',1, 'headerlines',2);
Ro = C{1}(1);

C = textscan(f_input,'%f',1, 'headerlines',2);
delta = C{1}(1);
Gamma = Ro * delta;


f_input = fopen([full_dir '/grid_def.all']);

linenum = 2;
C = textscan(f_input,'%s%s%u','delimiter',{' ','=','(',')'}, 'MultipleDelimsAsOne',1,'headerlines',linenum-1);
Nx = str2double(C{2}{1});
Ny = str2double(C{2}{2});
Nz = str2double(C{2}{3});
N_th = str2double(C{2}{4});

is3D = Ny > 1;


fclose(f_input);


if IC_type == 4 || IC_type == 6
    delta = 1;
    dTHdx = 1/(delta*Ro);
    dTHdy = 0; % In other horizontal direction
else
    dTHdx = 2/(Lx*Ro);
    dTHdy = 0;
end

dVdx = 0;
dVdz = dTHdx * delta*Ro;


% Read tau_crit

if exist([full_dir '/tau_crit'],'file')
    f_input = fopen([full_dir '/tau_crit']); % Units 1/f
    C = textscan(f_input,'%s',1);
    tau_crit = str2double(C{1}{1}) * Gamma; % GS Units
    fclose(f_input);
else
    tau_crit = 0;
end




