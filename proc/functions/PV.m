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



% [X,Z] = ndgrid(dx*[1:size(v,1)], cumsum(dzf));
% rdbu = flipud(cbrewer('div', 'RdBu', 128));
% figure()
% pcolor(X,Z,om_x); shading flat;
% colourbar('$\omega_x$',rdbu,'Centred');
% ylim([-0.01,0.04])
% xlabel('$x$')
% ylabel('$z$')
% 
% figure()
% pcolor(X,Z,om_z); shading flat;
% colourbar('$\omega_z$',rdbu,'Centred');
% ylim([-0.01,0.04])
% xlabel('$x$')
% ylabel('$z$')
% 
% figure()
% pcolor(X,Z,dxB); shading flat;
% colourbar('$\partial_x b$',rdbu,'Centred');
% ylim([-0.01,0.04])
% xlabel('$x$')
% ylabel('$z$')
% 
% figure()
% pcolor(X,Z,dzB); shading flat;
% colourbar('$\partial_z b$',rdbu,'Centred');
% ylim([-0.01,0.04])
% xlabel('$x$')
% ylabel('$z$')


end
