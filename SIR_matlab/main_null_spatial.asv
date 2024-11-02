% main_null_spatial.m
%                                                                                                                                    
% Package based off https://github.com/yingqiuz/SIR_simulator
% Zheng, Ying-Qiu, et al. PLoS biol. 17.11 (2019): e3000495.
%
% Gabriella Chan 12/10/23
% gabriella.chan@monash.edu
% Monash University
%
% We generate the simulated atrophy map for a gene pair, and compare with
% spatial nulls. Nulls are scrambles of empirical atrophy, simulated
% with BrainSMASH (brainsmash.readthedocs.io)
% We calculate the p-value of empirical atrophy against the null, and plot
% on a box-and-whiskers plot, with whiskers set to the 95% CI


function [] = main_null_spatial(clear, opt)
    % gene tables are indexed as T(:, ["gene1", "gene2"])
    arguments
        clear {mustBeText}                          % clear: clearance gene string array
        opt.risk {mustBeText}                       % risk: risk gene string array
        opt.parc {mustBeText} = "S132"              % accepts parcellations "DK", "S132", and "S332"
        opt.dbg {mustBeNumericOrLogical} = false;   % dbg: switches output to console or file
        opt.out {mustBeText} = "../SIR_results/results_3/gene_corrs.csv";
        opt.vis {mustBeNumericOrLogical} = false;   % vis: switches on visualisation options
        opt.null {mustBeText} = "none";
    end

    if opt.vis
        assert(~opt.null, 'This will produce too many figures.')
        asset(length(clear) * length(opt.risk) < 16, ...
            'This will produce too many figures.')
    end

    % load workspace variables in parameter object 
    params = SIRparameters(opt);
    params = set_atrophy(params);
    params = set_netw(params);
    ws_null = 'data/workspace_' + opt.parc + '_null.mat';
    params = set_null(params, opt, ws_null);
    % load gene expression data in genes parameter object
    genes = SIRgenes(opt, clear);

    % run SIR simulation for each gene pair
    genes = SIRiterator(params, genes);

    if (opt.dbg); display(genes.gene_corrs);
    else; write_async(genes.gene_corrs, opt.out); end
end

function [genes] = SIRiterator(params, genes)

    n = genes.n_risk * genes.n_clear;
    for idx = 1:n; genes = iterate(idx, params, genes); end

    % if params.vis
    %     % figure visualisation unavailable with parfor loops
    %     for idx = 1:n; genes = iterate(idx, params, genes); end
    % else
    %     % parfor idx = 1:n; gene_corrs(idx,:) = iterate(idx, params, genes); end
    %     parfor idx = 1:n; gene_corrs = [gene_corrs; iterate(idx, params, genes)]; end
    % end
    
    assert (~isempty(genes.gene_corrs));
    genes.gene_corrs = rmmissing(genes.gene_corrs);
end

% For each valid gene pair, simulate normal and misfolded protein growth 
% across all timepoints, and generate atrophy and correlation over time.
function genes = iterate (idx, params, genes)
    gene = SIRgene(genes, idx);
    if strcmp(gene.clear_name, gene.risk_name); return; end

    % SIRsimulator and SIRatrophy have visualization as a separate
    % parameter, switched off.
    [gene.Rnor_all, gene.Rmis_all] = SIRsimulator(params, gene, false);
    gene = SIRatrophy(params, gene, false);

    n = 1;
    if params.null == "spatial"; n = 1000; end
    for i = 1:n
        [bgs_max, cobre_max, hcpep_max, stages_max, tstep, ~] = ...
            SIRcorr(params, gene, i);
        avg_max = mean([bgs_max, cobre_max, hcpep_max]);
        
        row = (idx-1)*n + i;
        genes.gene_corrs(row, :) = {gene.risk_name, gene.clear_name, ...
            avg_max, bgs_max, cobre_max, hcpep_max, stages_max, tstep};
    end
end
