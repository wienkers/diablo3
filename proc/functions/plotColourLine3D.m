function f = plotColourLine3D(x,y,z,c, linewidth)
% Plots line in (x,y,z) with colours specified by c(x,y)

if nargin < 5
    linewidth = 2;
end

f = gcf;
% surface([x;x],[y;y],[z;z],[c;c],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',linewidth);

hold on;
for i = 1:size(x,2)
    patch([x(:,i); nan],[y(:,i); nan],[z(:,i); nan],[c(:,i); nan],'EdgeColor','interp','FaceColor','none','LineWidth',linewidth)
end

end