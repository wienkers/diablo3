function A = readField(filename, field, time_ind)

timename = sprintf('%04i',time_ind);

% Swap names for th in movies...
if ~strcmp(filename(end-6:end), 'mean.h5')
    if strcmp(field(1:2), 'b_')
        field = ['th1' field(2:end)];
    end

    if strcmp(field(1:2), 'bm')
        field = ['th' field(2:end)];
    end
end


field = ['/' field '/'];


A = h5read(filename,[field timename]);


end