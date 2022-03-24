function [kappa_t_rel, kappa_t_rel_mean, phi_d, phi_d_mean, Z_r, dTHdz_star]  = compute_KappaT(bpe_pdf,gradTH_pdf,   th_bin,  th_int,dTHdx,Lx,Lz, RePr)
%% Computes the turbulent Diffusivity relative to molecular (Kappa_t_rel), via Winters & D'Asaro 1996 Method
    
    Nbin = length(bpe_pdf);
    if size(bpe_pdf,1) ~= 1 && size(bpe_pdf,2) ~= 1
        error('BPE PDF should only be for 1 time!');
    end
    
    dTH = th_bin(2:end) - th_bin(1:end-1);
    
    dZ_add_to_bottom = th_int / (dTHdx*Lx);
    
    %% Smudge PDF in case the early-time is not random enough
    unfiltered_pdf = bpe_pdf;
%     for i = 2:Nbin-1
%         if unfiltered_pdf(i) == 0
%             bpe_pdf(i) = 0.25*(unfiltered_pdf(i+1) + unfiltered_pdf(i-1));
%             bpe_pdf(i+1) = bpe_pdf(i+1) - 0.25*unfiltered_pdf(i+1);
%             bpe_pdf(i-1) = bpe_pdf(i-1) - 0.25*unfiltered_pdf(i-1);
%         end
%     end
    
%     b = -3:3;
%     kern = exp(-b.^2/2^2)/sum(exp(-b.^2/2^2));
%     bpe_pdf = conv(bpe_pdf, kern, 'same');
    
%     if sum(bpe_pdf == 0) > 10
%         bpe_pdf = smooth(bpe_pdf,fix(Nbin/40));
%     end
    

    %% Compute Z_r of edges
    Z_r = zeros([Nbin+1,1]);
    Z_r(1) = dZ_add_to_bottom - 0.5*Lz;
    for i = 2:Nbin+1
        Z_r(i) = Z_r(i - 1) + bpe_pdf(i - 1) * Lz; % * dTH(i - 1)  !! PDF is the distribution NOT density -- so don't scale by dTH again...
    end
    
    
    %% Compute kappa_t
    % gradTH_pdf is a function of the th bins (which correspond now to Z_r)
    
    %dZ_r = Z_r(2:end) - Z_r(1:end-1);
    dZ_r = bpe_pdf;
    dTHdz_star = dTH ./ dZ_r;
    
%     mask = isinf(dTHdz_star) | isnan(dTHdz_star);
%     t    = 1:numel(dTHdz_star);
%     dTHdz_star(mask) = interp1(t(~mask), dTHdz_star(~mask), t(mask));
    
    
    kappa_t_rel = gradTH_pdf ./ dTHdz_star.^2;    
    
    phi_d = -1/RePr *  gradTH_pdf ./ dTHdz_star;
    
    
    %% Integrated Quantities
    
    
    phi_d_mean = sum(phi_d .* dZ_r); 
    kappa_t_rel_mean = sum(kappa_t_rel .* dZ_r);
    
    
    Z_r = Z_r(1:end-1).';


end

