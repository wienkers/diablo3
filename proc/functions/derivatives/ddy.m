function dfdy = ddy(f,dy)

N3 = size(f,3);

if N3 == 1 % x-y plane

    [Nx,Ny] = size(f);

    k2 = fftshift((-Ny/2 : Ny/2 - 1)*(2*pi/(Ny*dy)));
    ky = repmat(k2,[Nx,1]);
    fhat = fft(f,[],2);

    fhat = 1i*ky.*fhat;

    dfdy = ifft(fhat,[],2,'symmetric');
    
else
    
    [Nx,Nz,Ny] = size(f);

    k2 = fftshift((-Ny/2 : Ny/2 - 1)*(2*pi/(Ny*dy)));
    ky = repmat(reshape(k2,[1,1,Ny]),[Nx,Nz,1]);
    fhat = fft(f,[],3);

    fhat = 1i*ky.*fhat;

    dfdy = ifft(fhat,[],3,'symmetric');


end