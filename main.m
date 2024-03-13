% main.m
%                                                                                                                                   
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


function [] = main(clear_genes_in, risk_genes_in, debugging, plotting)
    
    % load gene expressions, real atrophy, ROIsize, functional connectivity...
    load('data/workspace_DK_pca.mat', 'genes', 'clear_score', 'risk_score', 'sconnDen', 'sconnLen', ...
        'ROIsize', 'emp_atrophy_bgs', 'emp_atrophy_cobre', 'emp_atrophy_hcpep', ...
        'emp_atrophy_stages');
    % load('data/workspace_DK_PD.mat', 'pd_atrophy_reordered', 'genes');

    all_clears = clear_score;
    all_risks = risk_score;

    % Input clearance genes and risk genes as string arrays
    
    % clearance and risk genes are input as gene x region expression tables
    % tables are indexed by variable names as T(:,["string", "array"])

    % if there are no input risk genes, all genes are taken as risk genes
    if ~exist("risk_genes_in", "var")
        risk_genes_in = string(all_risks.Properties.VariableNames);
    end
    if ~exist("debugging", "var")
        debugging = false;
        plotting = false;
    end

    clear_genes_bool = ismember(clear_genes_in, all_clears.Properties.VariableNames);
    valid_clear_genes = clear_genes_in(clear_genes_bool);
    clear_genes = all_clears(:, valid_clear_genes);

    risk_genes_bool = ismember(risk_genes_in, all_risks.Properties.VariableNames);
    valid_risk_genes = risk_genes_in(risk_genes_bool);
    risk_genes = all_risks(:,valid_risk_genes);

    if isempty(clear_genes) || isempty(risk_genes)
        disp(risk_genes_in)
        disp(all_risks(:, risk_genes_in))
        disp(clear_genes_in)
        disp(all_clears(:,clear_genes_in))
        disp('No valid gene combinations')
        return
    end
    
    N_regions = 41;
    v = 1;
    dt = 0.01;
    T_total = 20000;
    init_number = 1;
    prob_stay = 0.5;
    trans_rate = 1;
    % Single seed
    seeds = 1;
    seed = 40; % anterior hip
    % Multiple seeds
    % seeds = 41;
    index = 1;

    % Set up gene corr table
    % table is preallocated, so index for next row is passed through
    gene_corrs = table('Size', [width(risk_genes)*width(clear_genes)*seeds, 9], ...
    'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double'}, ...
    'VariableNames', {'risk gene', 'clearance gene', '~', 'correlation', 'bgs correlation', ...
    'cobre correlation', 'hcpep correlation', 'stages correlation', 'bgs t'});

    % Single seed
    [gene_corrs, ~] = SIRiterator(gene_corrs, index, 0, N_regions, v, dt, T_total, ...
        clear_genes, risk_genes, sconnLen, sconnDen, ROIsize, seed, ...
        init_number, prob_stay, trans_rate, emp_atrophy_bgs, ...
        emp_atrophy_cobre, emp_atrophy_hcpep, emp_atrophy_stages, plotting);

    % Multiple seeds
    % for seed = 1:seeds
    %     [gene_corrs, index] = SIRiterator(gene_corrs, index, N_regions, v, dt, T_total, ...
    %         clear_genes, risk_genes, sconnLen, sconnDen, ROIsize, seed, ...
    %         init_number, prob_stay, trans_rate, emp_atrophy_bgs, ...
    %         emp_atrophy_cobre, emp_atrophy_hcpep, emp_atrophy_stages, plotting);
    % end

    gene_corrs = rmmissing(gene_corrs);

    if (~isempty(gene_corrs))
        if (debugging)
            display(gene_corrs)
            % tail(gene_corrs)
        else
            out_file = 'results/gene_corrs_240311_pca.csv';
            write_async(gene_corrs, out_file)
        end
    end
end