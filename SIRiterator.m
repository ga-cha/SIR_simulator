% SIRiterator
%
% Gabriella Chan 28/06/23
% gabriella.chan@monash.edu
% Monash University
%
%
% This module serves as an intermediate layer to call SIR simulator
% and SIR atrophy.
%
% We iterate through all gene pairs and their gene expression and call SIR
% simulator to store an array of normal protein expression over regions 
% over time, and misfolded protein expression over regions over time.
%
% From this we derive theortical atrophy per region for each gene and
% correlate with empirical atrophy per region per gene.

function [gene_corrs, sim_atrophy] = SIRiterator(N_regions, v, dt, T_total, GBA, genes, ...
    sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, ...
    trans_rate, emp_atrophy)

gene_corrs = table('Size', [width(genes), 2], ...
    'VariableTypes', {'string', 'double'}, ...
    'VariableNames', {'gene', 'correlation'});

for gene = 1:width(genes)
    tic
    gene_name = genes.Properties.VariableNames{gene};
    disp(gene_name)
    %%%%% simulation ------ >>>
    [Rnor_all, Rmis_all] = SIRsimulator(N_regions, v, dt, T_total, GBA, ...
        genes(:,gene), sconnLen, sconnDen, ROIsize, seed, syn_control, ...
        init_number, prob_stay, trans_rate);
    [sim_atrophy] = SIRatrophy(Rnor_all, Rmis_all, sconnDen, N_regions, dt);
    [corr] = SIRcorr(sim_atrophy, emp_atrophy);
    gene_corrs(gene,:) = {gene_name, corr};
    toc
end

gene_corrs = sortrows(gene_corrs, 2);

end
