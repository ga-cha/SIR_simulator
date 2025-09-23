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

function sim_atrophy = sir_atrophy(params, proteins, vis)
    sc_weight = params.sc_weight;
    n_rois = params.n_rois;
    dt = params.dt;
    
    ratio = proteins.r_mis ./(proteins.r_nor + proteins.r_mis) ;
    ratio(proteins.r_mis < 1) = 0; % remove possible NaNs...
       
    % atrophy growth
    % GC: rescaled so deafferentation input is proportional to global 
    % region connectivity
    scaled_sc = sum(sc_weight)' ./ max(sum(sc_weight));
    % scaled_sc = 1;
    k2w = 0.5;
    k2 = k2w .* scaled_sc;
    k1 = 1 - k2w;

    % input weigths of deafferentation (scaled by structrual connectivity)
    weights = sc_weight ./ sum(sc_weight, 2);
    
    % neuronal loss caused by neighbouring deafferentation
    deaff = weights * (1-exp(-ratio * dt));
    % one time step back
    deaff = [zeros(n_rois, 1), deaff(:, 1:end-1)];
    ratio_cum = k1 .* (1-exp(-ratio * dt)) + k2 .* deaff;
    
    % add all the increments across t
    sim_atrophy = cumsum(ratio_cum, 2);

    % GC 
    sim_atrophy(sim_atrophy == 0) = eps;
    sim_atrophy = log(sim_atrophy);

    if vis && params.null == "none"
        plot_protein(proteins, sim_atrophy, params.roi_size);
    end

end

function plot_protein(proteins, sim_atrophy, roi_size)
    figure('Position', [150 150 1500 600], 'Color', 'white');
    subplot(1, 3, 1);
    plot((proteins.r_mis ./roi_size)');
    t = title ("Rmis/ROIsize dt");
    t.FontWeight = 'normal';
    xlabel("t")
    ylabel("atrophy")

    subplot(1, 3, 2);
    plot((proteins.r_mis ./(proteins.r_nor + proteins.r_mis))');
    t = title ("Rmis/Rnor+mis dt");
    t.FontWeight = 'normal';
    xlabel("t")
    ylabel("atrophy")

    subplot(1, 3, 3);
    plot((sim_atrophy)');
    t = title ("simulated atrophy dt");
    t.FontWeight = 'normal';
    xlabel("t")
    ylabel("atrophy")
end