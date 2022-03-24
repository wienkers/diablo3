function dfdz = ddz(f,dzf,BCs,varargin)
% For x-z plane. 2nd-order accurate
%   Additional arguments can define the:
%       Order of Derivative = 1 (default)

dzf = dzf(:)';

Ny = size(f,3);

if Ny ~=1 % 3D!
    dfdz = zeros(size(f));

    if nargin == 4
      order = fix(varargin{1});
    else
      order = 1;
    end

    for i = 1:Ny
        dfdz(:,:,i) = ddz(f(:,:,i),dzf,BCs,order);
    end
    return 
end


% Defaults
order = 1;
if nargin == 4
    order = fix(varargin{1});
end

if isnumeric(BCs)
    if BCs == 0
        BCs = 'Dirichlet';
    elseif BCs == 1
        BCs = 'Neumann';
    elseif BCs == 2
        BCs = 'Sided';
    else
        error('Not sure which BCs...');
    end
end


if (dzf(end) - dzf(1)) == 1
    error('You''re using gzf!')
end

if size(f,2) ~= length(dzf)
    flop = true;
    f = f.';
else
    flop = false;
end


if strcmp(BCs,'Sided') % Sided evaluation
    
    if order == 1
    
        dfdz = (f(:,3:end) - f(:,1:end-2)) ./ (2*dzf(2:end-1));
        dfdz = [(f(:,2)-f(:,1))/dzf(1) dfdz (f(:,end)-f(:,end-1))/dzf(end)];
        
    else
        error('Order not defined yet...');
    end
    
else

    if strcmp(BCs,'Dirichlet') % Dirichlet -- Odd extension

        fext = [-f(:,(order+1):-1:2)  f  -f(:,(end-1):-1:(end-order))];

    elseif strcmp(BCs,'Neumann') % Neumann -- Even extension

        fext = [f(:,(order+1):-1:2)  f  f(:,(end-1):-1:(end-order))];

    end
    
    if order == 1
        
        dfdz = (fext(:,3:end) - fext(:,1:end-2)) ./ (2*dzf);
        
    elseif order == 2
        
        fext = fext(:,2:end-1);
        dfdz = (fext(:,3:end) - 2*fext(:,2:end-1) + fext(:,1:end-2)) ./ dzf.^2;
        
    elseif order == 3
        
        fext = fext(:,2:end-1);
        dfdz = (0.5*fext(:,5:end) - fext(:,4:end-1) + fext(:,2:end-3) - 0.5*fext(:,1:end-4)) ./ dzf.^3;
        
    else
        error('Order not defined yet...');
    end
        
        
end




if flop
    dfdz = dfdz.';
end

end
