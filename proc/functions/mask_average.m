function f_avg = mask_average(f,mask)
% Averages just the masked quantity

f_avg = sum(f(mask)) / sum(mask(:));


end