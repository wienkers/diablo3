function dfdx = ddx(f,dx,varargin)
%% For x-z plane. 2nd-order accurate
%   Additional arguments can define either the:
%       Order of Derivative = 1 (default)
%       Derivative Type = 'Spectral' (default)


% Defaults
order = 1;
type = 'Spectral';

if nargin == 3
    if isnumeric(varargin{1})
        order = varargin{1};
    else
        type = varargin{1};
    end
elseif nargin == 4
    if isnumeric(varargin{1})
        order = varargin{1};
        type = varargin{2};
    else
        type = varargin{1};
        order = varargin{2};
    end
end



[Nx,Nz,Ny] = size(f);


if strcmp(type,'Spectral')
    
    if mod(Nx,2) == 1
        error('Nx must be even for Spectral Derivatives!')
    end
    
    k1 = fftshift((-Nx/2 : Nx/2 - 1)*(2*pi/(Nx*dx)));
    kx = repmat(k1',[1,Nz,Ny]);
    fhat = fft(f,[],1);

    fhat = (1i*kx).^order .* fhat;

    dfdx = ifft(fhat,[],1,'symmetric');
    
elseif strcmp(type,'FD')
    
    if order == 1
    
        dfdx = (f(3:end,:,:) - f(1:end-2,:,:)) / (2*dx);
        dfdx = [dfdx(1,:,:); dfdx; dfdx(end,:,:)];
        
    elseif order == 2
        
        dfdx = (f(3:end,:,:) - 2*f(2:end-1,:,:) + f(1:end-2,:,:)) / (dx^2);
        dfdx = [dfdx(1,:,:); dfdx; dfdx(end,:,:)];
        
    elseif order == 3
        
        dfdx = (0.5*f(5:end,:,:) - f(4:end-1,:,:) + f(2:end-3,:,:) - 0.5*f(1:end-4,:,:)) / (dx^3);
        dfdx = [dfdx(1,:,:); dfdx(1,:,:); dfdx; dfdx(end,:,:); dfdx(end,:,:)];
        
    else
        error('Derivative order not defined yet.')
    end 
    
else
    error('Derivative type not defined yet.');
end





end