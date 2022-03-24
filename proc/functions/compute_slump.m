function [xS_0, xS_sort] = compute_slump(A_b, x, Lx, Nz_avg)
% Computes the front slumped distance
% A_b must be full field

[~, Nz] = size(A_b);

%% Using straight averaging and 0-Crossing

btop = mean(A_b(:,1:Nz_avg),2);
bbot = mean(A_b(:,Nz-Nz_avg:end),2);

b = (btop - bbot(end:-1:1)) / 2;

xcross = b(1:end-2) < 0 & b(3:end) > 0;
xcross = [0; xcross; 0];

xS_0 = mean(x(xcross ~= 0)) - Lx/2;



%% Sort the top/bottom and then find the 0-Crossing

b_comp = [A_b(:,1:Nz_avg), -A_b(end:-1:1,Nz-Nz_avg:end)];

x_scaled = linspace(0,Lx,numel(b_comp)); x_scaled = x_scaled(1:end-1);
b_sort = sort(b_comp(:),'ascend');
xS_sort = x_scaled(b_sort(1:end-1) <  0 & b_sort(2:end) >  0) - Lx/2;



end
