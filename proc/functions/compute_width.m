function [width_5, width_tanh, b_sort_512] = compute_width(A_b, Lx)
% Computes the front width using sorting method
% A_b must be full field, scaled between -1 and 1

[Nx, Nz] = size(A_b);


%% Sort buoyancy (in x)
x_scaled = linspace(0,Lx,Nx*Nz);
b_sort = sort(A_b(:),'ascend')';


%% Width = x|b=0.45 - x|b=-0.45
x_scaled_m1 = x_scaled(1:end-1);
x_5n = x_scaled_m1(b_sort(1:end-1) < -0.45 & b_sort(2:end) > -0.45);
x_5p = x_scaled_m1(b_sort(1:end-1) <  0.45 & b_sort(2:end) >  0.45);

width_5 = (x_5p - x_5n); % Roughly corresponds to delta...



%% Width = Fit to tanh
fun = @(p,xdata) tanh((xdata - Lx/2) ./ p(1)); % Restrict the profile to be centred
width_tanh = lsqcurvefit(fun,1,x_scaled,b_sort);



%% Return also the sorted buoyancy but coarsened onto Nx = 512

x_s_512 = linspace(0,Lx,512);
b_sort_512 = interp1(x_scaled,b_sort,x_s_512);



end
