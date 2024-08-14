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

function [gene_corrs] = SIRiterator(params, genes)
    % Set up table to save gene correlation output
    gene_corrs = table(                                                     ...
        'Size', [width(genes.risk_genes) * width(genes.clear_genes), 8],    ...
        'VariableTypes', {'string', 'string', 'double', 'double', 'double', ...
        'double', 'double', 'double'},                                      ...
        'VariableNames', {'risk gene', 'clearance gene', 'correlation',     ...
        'bgs correlation', 'cobre correlation', 'hcpep correlation',        ...
        'stages correlation', 'bgs t'});

    n = genes.n_risk * genes.n_clear;
    
    if params.vis
        % switch to a for loop if visualising to the user
        for idx = 1:n; gene_corrs(idx,:) = iterate(idx, params, genes); end
    else
        parfor idx = 1:n; gene_corrs(idx,:) = iterate(idx, params, genes); end
    end
    

    % resize table for skipped rows
    gene_corrs = rmmissing(gene_corrs);
    assert (~isempty(gene_corrs));
end

% For each valid gene pair, simulate normal and misfolded protein growth 
% across all timepoints, and generate atrophy and correlation over time.
function [gene_corr] = iterate (idx, params, genes)
    % Rearrange index to  risk and clear gene indices
    [i, j] = ind2sub([genes.n_risk, genes.n_clear], idx);

    % the gene struct holds a risk/clear gene pair and intermediate outputs
    gene = struct();
    gene.risk_name = genes.risk_names.Value(i);
    gene.risk_gene = genes.risk_genes.Value(:, i);
    gene.clear_name = genes.clear_names.Value(j);
    gene.clear_gene = genes.clear_genes.Value(:, j);

    if strcmp(gene.clear_name, gene.risk_name); return; end

    % SIRsimulator and SIRatrophy have visualization as a separate
    % parameter, switched off.
    [gene.Rnor_all, gene.Rmis_all] = SIRsimulator(params, gene, false);
    gene = SIRatrophy(params, gene, false);
    [bgs_max, cobre_max, hcpep_max, stages_max, tstep, ~] = ...
        SIRcorr(params, gene);
    avg_max = mean([bgs_max, cobre_max, hcpep_max]);

    gene_corr = {gene.risk_name, gene.clear_name, ...
        avg_max, bgs_max, cobre_max, hcpep_max, stages_max, tstep};
end
