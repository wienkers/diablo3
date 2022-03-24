function [A_c,x_c,z_c] = coarsen(A,f,x,z)
%% Coarsens grid and data, by a factor f = [fx gz] in either direction

fx = f(1);
fz = f(2);

if nargin > 2
    x_c = conv(x,ones(fx,1)/fx,'valid');
        x_c = x_c(1:fx:end);
    z_c = conv(z,ones(fz,1)/fz,'valid');
        z_c = z_c(1:fz:end);
else
    x_c = [];
    z_c = [];
end


kernel = ones(fx, fz) / (fx*fz);

A = conv2(A,kernel,'valid');
    A_c = A(1:fx:end,1:fz:end);




end