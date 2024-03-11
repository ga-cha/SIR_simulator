function [] = main_null_rewired3(clear_genes_in, risk_genes_in, debuggy, plotting)
    % Input clearance genes and risk genes as string arrays
    if (~isstring(clear_genes_in)|~isstring(risk_genes_in))
        clear_genes_in = string(clear_genes_in);
        risk_genes_in = string(risk_genes_in);
    end

    % load gene expressions, real atrophy, ROIsize, functional connectivity...
    % simulated structural connectivity lengths are described in greater 
    % detail in Zheng PLoS biol. (2019)
    % The implementation is given in SIR_utils/null_ROI_dist.m
    load('data/workspace_rand100.mat', 'genes', 'sconnLen', ...
        'ROIsize', 'emp_atrophy_bgs', 'emp_atrophy_cobre', 'emp_atrophy_hcpep', ...
        'emp_atrophy_stages');
    % Load simulated structural connectivity density. Derived in
    % null_rewire.m, from BCT implementation of Maslov Sneppen
    load('data/rewire/sconnDen_sim.mat', 'sconnDen_sim');
    
    N_regions = 41;
    v = 1;
    dt = 0.01;
    T_total = 10000;
    init_number = 1;
    prob_stay = 0.5;
    trans_rate = 1;
    % initialise seed to hip
    seed = 40;


    % clearance and risk genes are input as gene x region expression tables
    % tables are indexed by variable names as T(:,["string", "array"])
    clear_genes_bool = ismember(clear_genes_in, genes.Properties.VariableNames);
    valid_clear_genes = clear_genes_in(clear_genes_bool);
    clear_genes = genes(:, valid_clear_genes);

    risk_genes_bool = ismember(risk_genes_in, genes.Properties.VariableNames);
    valid_risk_genes = risk_genes_in(risk_genes_bool);
    risk_genes = genes(:,valid_risk_genes);

    if isempty(clear_genes) || isempty(risk_genes)
        disp('No valid gene combinations')
        return
    end

    nulls = 10;
    % null_corrs = zeros(nulls,1);

    gene_corrs = table('Size', [width(risk_genes)*width(clear_genes)*nulls, 9], ...
    'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double'}, ...
    'VariableNames', {'risk gene', 'clearance gene', 'null', 'correlation', 'bgs correlation', ...
    'cobre correlation', 'hcpep correlation', 'stages correlation', 'bgs t'});
    index = 1;
    for i = 1:nulls
        [gene_corrs, index] = SIRiterator(gene_corrs, index, i, N_regions, v, dt, T_total, ...
            clear_genes, risk_genes, sconnLen, sconnDen_sim(:,:,i), ROIsize, seed, ...
            init_number, prob_stay, trans_rate, emp_atrophy_bgs, ...
        emp_atrophy_cobre, emp_atrophy_hcpep, emp_atrophy_stages, plotting);
        % null_corrs(i) = gene_corrs.correlation;
    end
    gene_corrs = sortrows(gene_corrs, 4);

    if (isempty(gene_corrs))
        if (debuggy)
            tail(gene_corrs)
        else
            out_file = 'results/gene_corrs_rand_rewired.csv';
            write_async(gene_corrs, out_file)
        end
    end
end