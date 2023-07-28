% Gabriella Chan 15/06/2023
% gabriella.chan@monash.edu
% Monash University

% A wrapper for calling SIR simulator
% iterates through a list of potential schizophrenia risk genes
% and clearance pathway markers (UPS or lysosomal degradation)

function [Rnor_all, Rmis_all, Rnor0, Pnor0, Pnor_all, Pmis_all] = SIRsimulator(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate)



%%%%% simulation ------ >>>
[Rnor_all, Rmis_all, Rnor0] = SIRsimulator(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate);
ratio = Rmis_all ./(Rnor_all + Rmis_all) ;
ratio(Rmis_all<1) = 0; % remove possible NaNs...

% atrophy growth
k1 = 0.5;
k2 = 1 - k1;
% input weigths of deafferentation (scaled by structrual connectivity)
weights = sconnDen ./ repmat(sum(sconnDen, 2), 1, N_regions);

% neuronal loss caused by lack of input from neighbouring regions
ratio_cum = weights * (1-exp(-ratio * dt));
% one time step back
ratio_cum = [zeros(N_regions, 1), ratio_cum(:, 1:end-1)];
ratio_cum = k2 * ratio_cum + k1 * (1-exp(-ratio * dt));

% add all the increments across t
simulated_atrophy = cumsum(ratio_cum, 2);
end
