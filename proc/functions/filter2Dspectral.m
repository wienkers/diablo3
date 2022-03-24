function A_f = filter2Dspectral(A,k_c, dx,gzf, BCs, symmetry)
%% Returns the filtered ("mean") field, exploiting x->-x, z->-z symmetry (where A->-A)
% BCs = 0 --> Dirichlet --> Sine
% BCs = 1 --> Neumann --> Cosine

[Nx,Nz] = size(A);
Lx = Nx*dx;
Lz = (gzf(end) - gzf(1));

Nx_c = 2 * ceil(k_c(1) / (2*pi/Lx)) + 2;
Nkz = Nz-2 + BCs; % For speed, want Nkz to be 2^m - 1 (for sine), or 2^m (for cosine)


% Use x-z symmetry to make average (mirror+flip)/2
if symmetry
    % Need to swap the left-most column of x to the right because the last
    %  point in x is no Lx, so this ensures the symmetry lines up!
    A = 0.5 * (A + -rot90([A(2:end,:); A(1,:)],2));
end


% Fourier transform each line across in x, but only keep first n modes

kx = fftshift((-Nx/2 : Nx/2 - 1) * (2*pi/Lx));
Fx = fft(A,[],1);
%Fx = Fx([1:Nx_c/2  end+1-Nx_c/2:end], :);
Fx((Nx_c/2+1):(Nx-Nx_c/2), :) = 0; % Only keep the modes < k_c

if BCs == 0 % Sine
    z_uni = (1 : Nkz) / (Nkz+1) * Lz;
    kz = (1 : Nkz) * pi/Lz;
else % Cosine
    z_uni = ((1 : Nkz) * 2 - 1) / (2*Nkz) * Lz;
    kz = (0 : Nkz-1) * pi/Lz;
end
dz_uni = z_uni(2) - z_uni(1);


for i = [1:Nx_c/2  (Nx+1-Nx_c/2):Nx]
    
    % Interpolate each Fourier coefficient onto uniform z grid
    Fx_interp = interp1(gzf,Fx(i,:),z_uni);
    
    % Cosine or Sine Transform in z depending on the BCs
    if BCs == 0 % Sine
        Fxz = dst(Fx_interp);
    else % Cosine
        Fxz = dct(Fx_interp);
    end
    
    % Set each of the coefficients for K outside of wavenumber ellipse
    K2 = (kx(i) / k_c(1))^2   +   (kz / k_c(2)).^2;
    Fxz(K2 > 1) = 0;
    
    % Inverse Cosine/Sine Transform
    if BCs == 0 % Sine
        Fx_interp = idst(Fxz);
    else % Cosine
        Fx_interp = idct(Fxz);
    end
    
    % Extend based on BCs & Interpolate back onto gzf
    if BCs == 0 % Sine/Odd
        % Boundary is exactly 0, then extend to negative
        z_uni_ext = [-dz_uni,   0,  z_uni,     Lz, (Lz+dz_uni)];
        Fx_interp_ext = [-Fx_interp(1), 0,  Fx_interp, 0,  -Fx_interp(end)];
    else % Cosine/Even
        % Boundary points are 0.5*dz away from either side
        z_uni_ext = [(-dz_uni+z_uni(1)),  z_uni,      (z_uni(end)+dz_uni)];
        Fx_interp_ext = [Fx_interp(1),          Fx_interp,  Fx_interp(end)];
    end
    
    Fx(i,:) = interp1(z_uni_ext, Fx_interp_ext, gzf);

end


% IFFT in x
A_f = ifft(Fx,[],1,'symmetric'); % Need to pad so that we get enough points in x!




end