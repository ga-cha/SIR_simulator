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

function gene = SIRatrophy(params, gene, vis)
    sconnDen = params.sconnDen;
    N_regions = params.N_regions;
    dt = params.dt;
    
    ratio = gene.Rmis_all ./(gene.Rnor_all + gene.Rmis_all) ;
    ratio(gene.Rmis_all<1) = 0; % remove possible NaNs...
       
    % atrophy growth
    k1 = 0.5;
    k2 = 1 - k1;
    % input weigths of deafferentation (scaled by structrual connectivity)
    weights = sconnDen ./ repmat(sum(sconnDen, 2), 1, N_regions);
    % % GC: An additional correction for (global) structural connectivity
    % weights = sconnDen .* rescale(sum(sconnDen)) ./ repmat(sum(sconnDen, 2), 1, N_regions);    
    
    % neuronal loss caused by lack of input from neighbouring regions
    ratio_cum = weights * (1-exp(-ratio * dt));
    % one time step back
    ratio_cum = [zeros(N_regions, 1), ratio_cum(:, 1:end-1)];
    ratio_cum = k2 * ratio_cum + k1 * (1-exp(-ratio * dt));
    
    % add all the increments across t
    gene.sim_atrophy = cumsum(ratio_cum, 2);

    if vis; plot_protein(gene.Rnor_all, gene.Rmis_all, params.ROIsize); end
end

function plot_protein(Rnor_all, Rmis_all, ROIsize)
    figure;
    plot((Rmis_all ./ROIsize)');
    t = title ("Rmis/ROIsize/t");
    t.FontWeight = 'normal';
    xlabel("t")
    ylabel("atrophy")

    figure;
    plot((Rmis_all ./(Rnor_all + Rmis_all))');
    t = title ("Rmis/Rnor+mis/t");
    t.FontWeight = 'normal';
    xlabel("t")
    ylabel("atrophy")

    figure;
    plot(((Rnor_all + Rmis_all)./ROIsize)');
    t = title ("Rnor+mis/ROIsize/t");
    t.FontWeight = 'normal';
    xlabel("t")
    ylabel("protein ratio")
end