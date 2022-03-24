function [A, xF, zF] = resample(A_n, x,z,f)
%% Resamples A with f as many points in both directions

if length(f) == 1
    f = [f f];
end

Nx = length(x);
Nz = length(z);

xF = linspace(min(x),max(x),f(1)*Nx);
zF = linspace(min(z),max(z),f(2)*Nz);
[XF,ZF] = ndgrid(xF,zF);

[X,Z] = ndgrid(x,z);

F = griddedInterpolant(X,Z,A_n.','spline');

A = F(XF,ZF).';



end