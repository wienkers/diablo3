function [A,time] = read3DField(filename, field)

time = h5readatt(filename,'/Timestep','Time');

if strcmp(field, 'b')
    field = 'TH1';
elseif strcmp(field, 'u')
    field = 'U';
elseif strcmp(field, 'v')
    field = 'V';
elseif strcmp(field, 'w')
    field = 'W'; % The name is swapped at run-time...
end


A = h5read(filename,['/Timestep/' field]);

end
