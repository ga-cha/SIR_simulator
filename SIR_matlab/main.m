% main.m
%
% Package based off https://github.com/yingqiuz/SIR_simulator
% Zheng, Ying-Qiu, et al. PLoS biol. 17.11 (2019): e3000495.
%
% Gabriella Chan 30/06/23
% gabriella.chan@monash.edu
% Monash University
%
% Simulates atrophy accrual due to the accumulation of misfolded
% protein aggregates, and correlates this with empirical data
%
% main("RNF144A", risk_names="RFXAP", dbg=true, vis=true);


function [] = main(clear_names, opt)
    % clr: list of clearance gene x region string arrays
    % rsk: list of risk gene x region string arrays
    % gene tables are indexed as T(:, ["gene1", "gene2"])
    % accepts parcellations "DK", "S132", and "S332", for DK + aseg,
    % Schaefer 100 + Tian S2 and Schaefer 300 + Tian S2 respectively
    % dbg: switches output to console or file
    % vis: switches on visualisation options
    arguments
        clear_names {mustBeText}
        opt.risk_names {mustBeText}
        opt.out {mustBeText} = '../../SIR_results/results_3/gene_corrs.csv';
        opt.parc {mustBeText} = "S132" ;
        opt.dbg {mustBeNumericOrLogical} = false;
        opt.vis {mustBeNumericOrLogical} = false;
        opt.null {mustBeText} = "none";
    end

    % load workspace variables into parameter object 
    params = SIR_parameters(opt);
    % load gene expression data into genes parameter object
    genes = SIR_genes(opt, clear_names);
    % start parallel pool
    p = gcp('nocreate'); if isempty(p); parpool('Threads'); end

    % tic
    gene_corrs = run_sim(genes, params);
    % toc

    if (opt.dbg)
        display(gene_corrs)
    else
        write_async(gene_corrs, opt.out)
    end
end


% We iterate through all gene pairs
% 
% 
% and their gene expression and call SIR
% simulator to store an array of normal protein expression over regions 
% over time, and misfolded protein expression over regions over time.
%
% From this we derive theortical atrophy per region for each gene and
% correlate with empirical atrophy per region per gene.

function gene_corrs = run_sim(genes, params)
    % run SIR simulation for each gene pair
    n = genes.n_risk * genes.n_clear;
    gene_objs = SIR_gene.empty(n, 0);
    
    if params.vis
        % figure visualisation unavailable with parfor loops
        for idx = 1:n
            gene = SIR_gene(genes, idx);
            gene_objs(idx) = gene.run_sim(params);
        end
    else
        parfor idx = 1:n
            gene = SIR_gene(genes, idx);
            gene_objs(idx) = gene.run_sim(params);
        end
    end

    gene_corrs = assemble_gene_corrs(params, gene_objs);
    assert (~isempty(gene_corrs));
end

% assemble table from gene pairs, preallocating for speed
function [gene_corrs] = assemble_gene_corrs(params, gene_objs)
    % collate correlation from each gene object into a single cell array
    gene_corr_cells = cell(length(gene_objs), 1);
    parfor i = 1:length(gene_objs)
        if isprop(gene_objs(i), 'gene_corr')
            gene_corr_cells{i} = gene_objs(i).gene_corr;
        end
    end
    % concatenate cell array into a single table
    C = vertcat(gene_corr_cells{:});
    gene_corrs = cell2table(C, 'VariableNames', [{'risk gene',               ...
        'clearance gene'}, params.emp_atr.Properties.VariableNames]);
end