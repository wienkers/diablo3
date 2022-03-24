function dfdt = ddt(f,time)
% For z-t plane

if size(f,2) ~= length(time)
    flop = true;
    f = f.';
else
    flop = false;
end

% Second Order...
dfdt = (f(:,3:end) - f(:,1:end-2)) ./ (time(3:end) - time(1:end-2));

dfdt = [  (f(:,2) - f(:,1)) ./ (time(2) - time(1))  ...
          dfdt ...
          (f(:,end) - f(:,end-1)) ./ (time(end) - time(end-1)) ];

% dfdt = (f(:,2:end) - f(:,1:end-1)) ./ (time(2:end) - time(1:end-1));
% dfdt = [dfdt dfdt(:,end)];



if flop
    dfdt = dfdz.';
end


end