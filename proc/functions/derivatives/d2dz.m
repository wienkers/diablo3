function d2fdz = d2dz(f,dzf,BCs)
% For x-z plane
%% BEING DEPRECATED

if (dzf(end) - dzf(1)) == 1
    error('You''re using gzf!')
end

if size(f,2) ~= length(dzf)
    flop = true;
    f = f.';
else
    flop = false;
end

if BCs == 0 % Dirichlet -- Odd extension
    fext = [-f(:,2)  f  -f(:,end-1)];
elseif BCs == 1 % Neumann -- Even extension
    fext = [f(:,2)  f  f(:,end-1)];
end

d2fdz = (fext(:,3:end) - 2*fext(:,2:end-1) + fext(:,1:end-2)) ./ dzf.^2;


if flop
    d2fdz = d2fdz.';
end


end