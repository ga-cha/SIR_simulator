% SIRcorr
%
% Gabriella Chan 28/06/23
% gabriella.chan@monash.edu
% Monash University
%
% Given simulated atrophy, this module calculates similarity with
% empirical atrophy.
% Simulated atrophy input contains atrophy per region per timepoint. For
% each timepoint, then, we can calculate empirical correlation. We then
% return the maximum.

function [corr, tstep] = SIRcorr(sim_atrophy, emp_atrophy)
    % something like: for each timepoint, pair atrophy values. 
    %for roi = 1:length(sim_atrophy)
    % Calculate the correlation coefficient
    max_corr = zeros(2, 2);
    tstep = 0;
    for i = 1:20000
        corr_mat = corrcoef(sim_atrophy(:, i), emp_atrophy);
        if max_corr(1, 2) < corr_mat (1, 2)
            max_corr = corr_mat;
            tstep = i;
        end
    end
    % extract the correlation coefficient
    corr = max_corr(1, 2);
    %end
end