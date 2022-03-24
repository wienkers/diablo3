function A = filter2D(A_n, r, dx, dzf)
%% Filters data in A_n on the X,Z grid using Gaussian radius r

if length(r) == 1
    r = [r r];
end

wN_x = fix(6*r(1)/dx);
if mod(wN_x,2) == 0
    wN_x = wN_x + 1;
end

dz = mean(dzf);

wN_z = fix(6*r(2)/dz);
if mod(wN_z,2) == 0
    wN_z = wN_z + 1;
end

x = -dx*(wN_x-1)/2:dx:dx*(wN_x-1)/2;
z = -dz*(wN_z-1)/2:dz:dz*(wN_z-1)/2;
[X, Z] = ndgrid(x,z);


% Gaussian kernel
kernel = exp(-(X.^2/r(1)^2 + Z.^2/r(2)^2));
kernel = kernel / sum(kernel(:)); %trapz(x,trapz(z,kernel,2),1);

N = size(kernel);
A_n_pad = padarray(A_n,N,'replicate','both');

A = conv2(A_n_pad,kernel,'same');

A = A(N(1):end-N(1)-1,N(2):end-N(2)-1);

end