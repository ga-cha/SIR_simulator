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

function [corr] = SIRcorr(sim_atrophy, emp_atrophy)
    % something like: for each timepoint, pair atrophy values. 
    %for roi = 1:length(sim_atrophy)
    % Calculate the correlation coefficient
    corr_mat = corrcoef(sim_atrophy(:, 20000), emp_atrophy);
    % extract the correlation coefficient
    corr = corr_mat(1, 2);
    %end
end