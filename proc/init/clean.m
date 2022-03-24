% Clears all variables except for those read by read_mean.m

if exist('mean_variables','var')
    clearvars('-except',mean_variables{:});
end
