% SIRatrophy
%
% Gabriella Chan 28/06/23
% gabriella.chan@monash.edu
% Monash University
%
% We take normal and misfolded protein expression over time and return
% simulated atrophy. Atrophy may occur as the result of misfolded protein
% accumulation, or deafferentiation.
% Simulated atrophy is returned as atrophy per region per timepoint
%

function [sim_atrophy] = SIRatrophy(Rnor_all, Rmis_all, ROIsize, sconnDen, N_regions, dt, plotting)

ratio = Rmis_all ./(Rnor_all + Rmis_all) ;
% ratio = Rmis_all ./ROIsize;
ratio(Rmis_all<1) = 0; % remove possible NaNs...

if plotting
    % figure;
    % plot((Rmis_all ./ROIsize)');
    % t = title ("Rmis/ROIsize/t");
    % t.FontWeight = 'normal';
    % xlabel("t")
    % ylabel("atrophy")
    % 
    % figure;
    % plot((Rmis_all ./(Rnor_all + Rmis_all))');
    % t = title ("Rmis/Rnor+mis/t");
    % t.FontWeight = 'normal';
    % xlabel("t")
    % ylabel("atrophy")
    % 
    % figure;
    % plot(((Rnor_all + Rmis_all)./ROIsize)');
    % t = title ("Rnor+mis/ROIsize/t");
    % t.FontWeight = 'normal';
    % xlabel("t")
    % ylabel("protein ratio")
end

% atrophy growth
k1 = 0.5;
k2 = 1 - k1;
% input weigths of deafferentation (scaled by structrual connectivity)
weights = sconnDen ./ repmat(sum(sconnDen, 2), 1, N_regions);

% % GC: An additional correction for (global) structural connectivity (and
% region size)
% weights = sconnDen .* rescale(sum(sconnDen)./ROIsize') ./ repmat(sum(sconnDen, 2), 1, N_regions);
% weights = sconnDen .* rescale(sum(sconnDen)) ./ repmat(sum(sconnDen, 2), 1, N_regions);


% neuronal loss caused by lack of input from neighbouring regions
ratio_cum = weights * (1-exp(-ratio * dt));
% one time step back
ratio_cum = [zeros(N_regions, 1), ratio_cum(:, 1:end-1)];
ratio_cum = k2 * ratio_cum + k1 * (1-exp(-ratio * dt));

% add all the increments across t
sim_atrophy = cumsum(ratio_cum, 2);

end