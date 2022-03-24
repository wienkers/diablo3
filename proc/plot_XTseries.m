%% Plots z = 0 TimeSeries of various quantities, now for *finite* fronts
%   NOTE: Uses Mean Slices Now
close all; clean;

set_inputs;
make_grids;
init_movie;

mkdir(full_dir, '/timeseries');

bmidz = zeros([Nt,Nx]);
umidz = zeros([Nt,Nx]);
vmidz = zeros([Nt,Nx]);
wmidz = zeros([Nt,Nx]);

bbotz = zeros([Nt,Nx]);
ubotz = zeros([Nt,Nx]);
vbotz = zeros([Nt,Nx]);
wbotz = zeros([Nt,Nx]);

btopz = zeros([Nt,Nx]);  
utopz = zeros([Nt,Nx]);
vtopz = zeros([Nt,Nx]);
wtopz = zeros([Nt,Nx]);

b8topz = zeros([Nt,Nx]);
u8topz = zeros([Nt,Nx]);
v8topz = zeros([Nt,Nx]);
w8topz = zeros([Nt,Nx]);

b8botz = zeros([Nt,Nx]);
u8botz = zeros([Nt,Nx]);
v8botz = zeros([Nt,Nx]);
w8botz = zeros([Nt,Nx]);

t = zeros([1,Nt]);


for k = 1:Nt
    timename = sprintf('%04i',k);

    t(k) = h5readatt(filename,['/u_xz/' timename],'Time');

    bmidz(k,:) = h5read(filename_mean,['/thme_xz/' timename],[1,fix(Nz/2)],[Nx,1],[1,1]);
    b8botz(k,:) = h5read(filename_mean,['/thme_xz/' timename],[1,fix(Nz/8)],[Nx,1],[1,1]);
    bbotz(k,:) = h5read(filename_mean,['/thme_xz/' timename],[1,4],[Nx,1],[1,1]);
    b8topz(k,:) = h5read(filename_mean,['/thme_xz/' timename],[1,Nz-fix(Nz/8)],[Nx,1],[1,1]);
    btopz(k,:) = h5read(filename_mean,['/thme_xz/' timename],[1,Nz-4],[Nx,1],[1,1]);

    umidz(k,:) = h5read(filename_mean,['/ume_xz/' timename],[1,fix(Nz/2)],[Nx,1],[1,1]);
    u8botz(k,:) = h5read(filename_mean,['/ume_xz/' timename],[1,fix(Nz/8)],[Nx,1],[1,1]);
    ubotz(k,:) = h5read(filename_mean,['/ume_xz/' timename],[1,4],[Nx,1],[1,1]);
    u8topz(k,:) = h5read(filename_mean,['/ume_xz/' timename],[1,Nz-fix(Nz/8)],[Nx,1],[1,1]);
    utopz(k,:) = h5read(filename_mean,['/ume_xz/' timename],[1,Nz-4],[Nx,1],[1,1]);

    vmidz(k,:) = h5read(filename_mean,['/vme_xz/' timename],[1,fix(Nz/2)],[Nx,1],[1,1]);
    v8botz(k,:) = h5read(filename_mean,['/vme_xz/' timename],[1,fix(Nz/8)],[Nx,1],[1,1]);
    vbotz(k,:) = h5read(filename_mean,['/vme_xz/' timename],[1,4],[Nx,1],[1,1]);
    v8topz(k,:) = h5read(filename_mean,['/vme_xz/' timename],[1,Nz-fix(Nz/8)],[Nx,1],[1,1]);
    vtopz(k,:) = h5read(filename_mean,['/vme_xz/' timename],[1,Nz-4],[Nx,1],[1,1]);
    
    wmidz(k,:) = h5read(filename_mean,['/wme_xz/' timename],[1,fix(Nz/2)],[Nx,1],[1,1]);
    w8botz(k,:) = h5read(filename_mean,['/wme_xz/' timename],[1,fix(Nz/8)],[Nx,1],[1,1]);
    wbotz(k,:) = h5read(filename_mean,['/wme_xz/' timename],[1,4],[Nx,1],[1,1]);
    w8topz(k,:) = h5read(filename_mean,['/wme_xz/' timename],[1,Nz-fix(Nz/8)],[Nx,1],[1,1]);
    wtopz(k,:) = h5read(filename_mean,['/wme_xz/' timename],[1,Nz-4],[Nx,1],[1,1]);
    
    


end


dxbmidz = ddx(bmidz.',dx).' + dTHdx;
dxb8botz = ddx(b8botz.',dx).' + dTHdx;
dxbbotz = ddx(bbotz.',dx).' + dTHdx;
dxb8topz = ddx(b8topz.',dx).' + dTHdx;
dxbtopz = ddx(btopz.',dx).' + dTHdx;


dxvmidz = ddx(vmidz.',dx).' + dVdx;
dxv8botz = ddx(v8botz.',dx).' + dVdx;
dxvbotz = ddx(vbotz.',dx).' + dVdx;
dxv8topz = ddx(v8topz.',dx).' + dVdx;
dxvtopz = ddx(vtopz.',dx).' + dVdx;


dxumidz = ddx(umidz.',dx).';
dxu8botz = ddx(u8botz.',dx).';
dxubotz = ddx(ubotz.',dx).';
dxu8topz = ddx(u8topz.',dx).';
dxutopz = ddx(utopz.',dx).';



% Derived & Smoothed
N_smooth = 10;

Gammamidz = smoothdata(Gamma*dxbmidz, 2,'gaussian',N_smooth);
Gamma8botz = smoothdata(Gamma*dxb8botz, 2,'gaussian',N_smooth);
Gammabotz = smoothdata(Gamma*dxbbotz, 2,'gaussian',N_smooth);
Gamma8topz = smoothdata(Gamma*dxb8topz, 2,'gaussian',N_smooth);
Gammatopz = smoothdata(Gamma*dxbtopz, 2,'gaussian',N_smooth);

f_effmidz = real(sqrt(1 + smoothdata(dxvmidz * Gamma, 2,'gaussian',N_smooth) ));
f_eff8botz = real(sqrt(1 + smoothdata(dxv8botz * Gamma, 2,'gaussian',N_smooth) ));
f_effbotz = real(sqrt(1 + smoothdata(dxvbotz * Gamma, 2,'gaussian',N_smooth) ));
f_eff8topz = real(sqrt(1 + smoothdata(dxv8topz * Gamma, 2,'gaussian',N_smooth) ));
f_efftopz = real(sqrt(1 + smoothdata(dxvtopz * Gamma, 2,'gaussian',N_smooth) ));

zetamidz = smoothdata(dxvmidz * Gamma, 2,'gaussian',N_smooth);
zeta8botz = smoothdata(dxv8botz * Gamma, 2,'gaussian',N_smooth);
zetabotz = smoothdata(dxvbotz * Gamma, 2,'gaussian',N_smooth);
zeta8topz = smoothdata(dxv8topz * Gamma, 2,'gaussian',N_smooth);
zetatopz = smoothdata(dxvtopz * Gamma, 2,'gaussian',N_smooth);

h_convmidz = smoothdata(dxumidz * Gamma, 2,'gaussian',N_smooth);
h_conv8botz = smoothdata(dxu8botz * Gamma, 2,'gaussian',N_smooth);
h_convbotz = smoothdata(dxubotz * Gamma, 2,'gaussian',N_smooth);
h_conv8topz = smoothdata(dxu8topz * Gamma, 2,'gaussian',N_smooth);
h_convtopz = smoothdata(dxutopz * Gamma, 2,'gaussian',N_smooth);

wmidz = smoothdata(wmidz, 2,'gaussian',N_smooth);
w8botz = smoothdata(w8botz, 2,'gaussian',N_smooth);
wbotz = smoothdata(wbotz, 2,'gaussian',N_smooth);
w8topz = smoothdata(w8topz, 2,'gaussian',N_smooth);
wtopz = smoothdata(wtopz, 2,'gaussian',N_smooth);


x = x';

tl = max(t)/Gamma;

if IC_type == 5 || IC_type == 7
  % Make full field for b
  for i = 1:size(bmidz,1)
      bmidz(i,:) = bmidz(i,:) + dTHdx*(x-Lx/2); % To get full buoyancy field
      b8botz(i,:) = b8botz(i,:) + dTHdx*(x-Lx/2);
      bbotz(i,:) = bbotz(i,:) + dTHdx*(x-Lx/2);
      b8topz(i,:) = b8topz(i,:) + dTHdx*(x-Lx/2);
      btopz(i,:) = btopz(i,:) + dTHdx*(x-Lx/2);
  end

  % Make full velocity field for v
  vmidz = vmidz + dVdz * (gzf(fix(Nz/2)));% - Lz/2);
  v8botz = v8botz + dVdz * (gzf(fix(Nz/8)));% - Lz/2);
  vbotz = vbotz + dVdz * (gzf(4));% - Lz/2);
  v8topz = v8topz + dVdz * (gzf(Nz-fix(Nz/8)));% - Lz/2);
  vtopz = vtopz + dVdz * (gzf(Nz-4));% - Lz/2);

end






f = plotTimeSeries(t/Gamma,x,bmidz.','$b(z=1/2)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/b_midZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,Gammamidz.','$\Gamma(z=1/2) / \Gamma_0$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/dxb_midZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,umidz.','$u(z=1/2)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/u_midZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,vmidz.','$v(z=1/2)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/v_midZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,wmidz.','$\tilde{w}(z=1/2)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/w_midZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,f_effmidz.','');  xlim([0,tl]);
colourbar('$\tilde{f}_\mathrm{eff}(z=1/2)/f$',gn,'Percentile98Positive');
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/feff_midZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,zetamidz.','$\zeta(z=1/2)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/zeta_midZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,h_convmidz.','$\nabla_h u(z=1/2)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/hconv_midZ_XT.png'])



f = plotTimeSeries(t/Gamma,x,b8botz.','$b(z=1/8)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/b_midbotZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,Gamma8botz.','$\Gamma(z=1/8) / \Gamma_0$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/dxb_midbotZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,u8botz.','$u(z=1/8)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/u_midbotZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,v8botz.','$v(z=1/8)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/v_midbotZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,w8botz.','$\tilde{w}(z=1/8)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/w_midbotZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,f_eff8botz.','');  xlim([0,tl]);
colourbar('$\tilde{f}_\mathrm{eff}(z=1/8)/f$',gn,'Percentile98Positive');
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/feff_midbotZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,zeta8botz.','$\zeta(z=1/8)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/zeta_midbotZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,h_conv8botz.','$\nabla_h u(z=1/8)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/hconv_midbotZ_XT.png'])



f = plotTimeSeries(t/Gamma,x,bbotz.','$b(z=4\Delta z)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/b_botZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,Gammabotz.','$\Gamma(z=4\Delta z) / \Gamma_0$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/dxb_botZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,ubotz.','$u(z=4\Delta z)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/u_botZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,vbotz.','$v(z=4\Delta z)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/v_botZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,wbotz.','$\tilde{w}(z=4\Delta z)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/w_botZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,f_effbotz.','');  xlim([0,tl]);
colourbar('$\tilde{f}_\mathrm{eff}(z=4\Delta z)/f$',gn,'Percentile98Positive');
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/feff_botZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,zetabotz.','$\zeta(z=4\Delta z)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/zeta_botZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,h_convbotz.','$\nabla_h u(z=4\Delta z)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/hconv_botZ_XT.png'])




f = plotTimeSeries(t/Gamma,x,b8topz.','$b(z=7/8)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/b_midtopZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,Gamma8topz.','$\Gamma(z=7/8) / \Gamma_0$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/dxb_midtopZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,u8topz.','$u(z=7/8)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/u_midtopZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,v8topz.','$v(z=7/8)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/v_midtopZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,w8topz.','$\tilde{w}(z=7/8)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/w_midtopZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,f_eff8topz.','');  xlim([0,tl]);
colourbar('$\tilde{f}_\mathrm{eff}(z=7/8)/f$',gn,'Percentile98Positive');
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/feff_midtopZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,zeta8topz.','$\zeta(z=7/8)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/zeta_midtopZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,h_conv8topz.','$\nabla_h u(z=7/8)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/hconv_midtopZ_XT.png'])



f = plotTimeSeries(t/Gamma,x,btopz.','$b(z=H-4\Delta z)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/b_topZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,Gammatopz.','$\Gamma(z=H-4\Delta z) / \Gamma_0$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/dxb_topZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,utopz.','$u(z=H-4\Delta z)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/u_topZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,vtopz.','$v(z=H-4\Delta z)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/v_topZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,wtopz.','$\tilde{w}(z=H-4\Delta z)$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/w_topZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,f_efftopz.','');  xlim([0,tl]);
colourbar('$\tilde{f}_\mathrm{eff}(z=H-4\Delta z)/f$',gn,'Percentile98Positive');
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/feff_topZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,zetatopz.','$\zeta(z=H-4\Delta z)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/zeta_topZ_XT.png'])

f = plotTimeSeries(t/Gamma,x,h_convtopz.','$\nabla_h u(z=H-4\Delta z)/f$');  xlim([0,tl]);
ylabel('$x$')
save_fig(f,[full_dir '/timeseries/hconv_topZ_XT.png'])



clear f;
save([full_dir '/XT_series.mat']);
