function t = read_time(filename,k)

if strcmp(filename(end-6:end), 'mean.h5')
    t = h5read(filename,['/time/' sprintf('%04i',k)]);
else % Movie file...
    t = h5readatt(filename,['/u_xz/' sprintf('%04i',k)],'Time');
end

end