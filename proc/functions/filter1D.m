function A = filter1D(A_n, r, dx)
%% Filters data in A_n with Gaussian radius r

dx = mean(dx);

wN_x = fix(6*r/dx);
if mod(wN_x,2) == 0
    wN_x = wN_x + 1;
end

x = -dx*(wN_x-1)/2:dx:dx*(wN_x-1)/2;

% Gaussian kernel
kernel = exp(-(x.^2) / (r^2));
kernel = kernel / sum(kernel(:));


A = conv(A_n,kernel,'same');

end