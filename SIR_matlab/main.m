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
        opt.parc {mustBeText} = "S132" 
        opt.dbg {mustBeNumericOrLogical} = false;
        opt.vis {mustBeNumericOrLogical} = false;
        opt.null {mustBeText} = "none";
    end

    % load workspace variables into parameter object 
    params = SIR_parameters(opt);
    params = set_atrophy(params);
    params = set_netw(params);
    % load gene expression data into genes parameter object
    genes = SIR_genes(opt, clear_names);
    % start parallel pool
    p = gcp('nocreate'); if isempty(p); parpool('Threads'); end

    tic
    % run SIR simulation for each gene pair
    gene_corrs = sir_iterator(params, genes);
    toc

    if (opt.dbg)
        display(gene_corrs)
    else
        write_async(gene_corrs, opt.out)
    end
end

