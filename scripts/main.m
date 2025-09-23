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

function gene_corrs = main(clear_names, opt)
    % clear_names: clearance gene name string array
    % risk_names:  risk gene name string array
    % parc: accepts "DK", "S132", and "S332" only, for DK + aseg,
    %       Schaefer 100 + Tian S2 and Schaefer 300 + Tian S2
    % seed: seed region for misfolded protein injection
    % dbg:  debugging flag, switches output to console or file
    % vis:  visualisation flag
    % pf:   parfor flag
    % null: accepts "none", "spatial", and "rewired"
    arguments
        clear_names {mustBeText}
        opt.risk_names {mustBeText}
        opt.out {mustBeText} = '../results/gene_corrs.csv';
        opt.parc {mustBeText} = "S132" ;
        opt.seed {mustBeNonnegative}
        opt.init {mustBeNonnegative} = 1;
        opt.dbg {mustBeNumericOrLogical} = false;
        opt.vis {mustBeNumericOrLogical} = false;
        opt.pf {mustBeNumericOrLogical} = true;
        opt.null {mustBeText} = "none";
    end

    % load workspace variables into parameter object 
    params = SIR_parameters(opt);
    % load gene expression data into genes parameter object
    genes = SIR_genes(opt, clear_names);
    % start parallel pool
    p = gcp('nocreate'); 
    if (isempty(p) && params.pf == true)
        parpool('Threads'); 
    end

    % tic
    gene_corrs = run_sim(genes, params);
    % toc

    if (opt.dbg)
        display(gene_corrs)
    else
        write_async(gene_corrs, opt.out)
    end
end


% Iterates through all gene pairs and generates gene pair object.
% Gene correlations are assembled from collated list of gene objects
function gene_corrs = run_sim(genes, params)
    % run SIR simulation for each gene pair
    n = genes.n_risk * genes.n_clear;
    gene_objs = SIR_gene.empty(n, 0);
    if params.pf
        parfor idx = 1:n
            gene = SIR_gene(genes, params, idx);
            gene_objs(idx) = gene.run_gene(params);
        end
    else
        for idx = 1:n
            gene = SIR_gene(genes, params, idx);
            gene_objs(idx) = gene.run_gene(params);
        end
    end

    gene_corrs = assemble_gene_corrs(params, gene_objs);
    assert (~isempty(gene_corrs));
end

% assemble table from gene pairs
function [gene_corrs] = assemble_gene_corrs(params, gene_objs)
    % collate correlation from each gene object into a single cell array
    gene_corr_cells = cell(length(gene_objs), 1);
    for i = 1:length(gene_objs)
        if isprop(gene_objs(i), 'gene_corr')
            gene_corr_cells{i} = gene_objs(i).gene_corr;
        end
    end
    % concatenate cell array into a single table
    C = vertcat(gene_corr_cells{:});
    if params.null == "none"
        gene_corrs = cell2table(C, "VariableNames", [gene_objs(1).col_names, ...
            params.emp_atr.Properties.VariableNames]);
    else
        gene_corrs = cell2table(C, "VariableNames", [gene_objs(1).col_names, ...
            params.emp_atr.Properties.VariableNames, "p_" + params.null]);
        gene_corrs = removevars(gene_corrs, ["t_max", "t_05"]);
    end
end