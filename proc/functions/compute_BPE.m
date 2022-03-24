function [BPE_full, Z_r, dbdz_s]  = compute_BPE(bpe_pdf, th_bin,  b_cutoff_high, b_cutoff_low,   th_int,dTHdx,Lx,Lz)
%% Computes the BPE from the sorted PDF
%  Can also do reduced BPE (i.e. between isopycnals.
    
    Nbin = length(bpe_pdf);
    if size(bpe_pdf,1) ~= 1 && size(bpe_pdf,2) ~= 1
        error('BPE PDF should only be for 1 time!');
    end
    
    dTH = th_bin(2:end) - th_bin(1:end-1);
    
    th_bin_centre = 0.5*(th_bin(2:end) + th_bin(1:end-1));
    
    
    % Determine the cutoff indices
    [~,ind_start_pdf] = min(abs(th_bin_centre - b_cutoff_low));
    [~,ind_end_pdf] = min(abs(th_bin_centre - b_cutoff_high));
    
    
    %% Re-scale PDF (since it was scaled as a density function at run-time...)
    %bpe_pdf = bpe_pdf / sum(bpe_pdf); % NOTE: This must be done in read_mean...

    
    dZ_add_to_bottom = th_int / (dTHdx*Lx);
    
    %% Compute Z_r of edges
    Z_r = zeros([1,Nbin+1]);
    Z_r(1) = dZ_add_to_bottom - 0.5*Lz;
    for i = 2:Nbin+1
        Z_r(i) = Z_r(i - 1) + bpe_pdf(i - 1) * Lz; % * dTH(i - 1)  !! PDF is the distribution NOT density -- so don't scale by dTH again...
    end    
    
    
    Z_c = 0.5 * (Z_r(1:end-1) + Z_r(2:end)); % Centred Z_r !
    
    
    Z_r_upper = Z_r(ind_end_pdf+1);
    Z_r_lower = Z_r(ind_start_pdf);
    Lz_reduced = Z_r_upper - Z_r_lower;
    
        
    %% Compute BPE
    BPE_full = 0;
    for i = 1:Nbin  % Integrate
        BPE_full = BPE_full - (th_bin_centre(i) * Z_c(i)) * (Z_r(i + 1) - Z_r(i)) / Lz;
    end
    
    
    if b_cutoff_high ~= 0 && b_cutoff_low ~= 0
        BPE_iso = 0;
        for i = ind_start_pdf:ind_end_pdf  % Integrate
            BPE_iso = BPE_iso - (th_bin_centre(i) * Z_c(i)) * (Z_r(i + 1) - Z_r(i)) / Lz; %Lz_reduced;
        end
    end
    
    
    [~,ind_middle] = min(abs(Z_c));
    
    dbdz_s = (th_bin(ind_middle+1) - th_bin(ind_middle)) / (Z_r(ind_middle+1) - Z_r(ind_middle));
    
    
    
    
%     %% BPE along slope of Gamma
% 
%     Z_v = Z_r * Gamma / sqrt(1+Gamma^2);
% 
%     SI_BPE(k) = 0;
%     for i = ind_start_pdf:ind_end_pdf  % Integrate
%         SI_BPE(k) = SI_BPE(k) - (th_bins_centre(i) * 0.5 * (Z_v(i + 1) + Z_v(i))) * (Z_r(i + 1) - Z_r(i)) / Lz_reduced;
%     end



end

