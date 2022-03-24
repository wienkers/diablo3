function [q, q_BC, q_vort] = PV(v,b,RoDelta,dTHdx,dVdx,dVdz,dx,dzf,gzf)
% For x-z plane
% Inputs are _perturbation_ quantities as output from Diablo
%  Here we assume that BCs are stress free and insulating
%    so that ddz(v_full) and ddz(b_full) are 0

gzf = gzf(:)';
Z = repmat(gzf,[size(v,1), 1]);
v_full = v + dVdz*Z;

om_x = -ddz(v_full,dzf,1);%  - dVdz;  % *Need* to use v_full to have even extension!
%om_y = ddz(u,dzf,1) - ddx(w,dx);
om_z = ddx(v,dx) + dVdx;

dxB = ddx(b,dx) + dTHdx;

dzB = ddz(b,dzf,1);

q_BC = om_x.*dxB;
q_vort = (om_z + 1/(RoDelta)).*dzB;

q = q_BC + q_vort; %om_x.*dxB + (om_z + 1/(RoDelta)).*dzB;


end
