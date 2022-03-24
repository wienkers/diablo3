function f = plotColourLine(x,y,z,linewidth)
% Plots y vs x but colours the line based on z(x,y)

if nargin < 4
    linewidth = 2;
end

f = gcf;
surface([x;x],[y;y],[zeros(size(x));zeros(size(x))],[z;z],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',linewidth);


end