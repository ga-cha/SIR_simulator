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

function [] = main2(clear_genes_in, risk_genes_in, plotting)
    assert(isstring(clear_genes_in) && isstring(risk_genes_in), ...
        "Input clearance genes and risk genes as string arrays")
    % load gene expressions, real atrophy, ROIsize, functional connectivity...
    load('data/workspace_DK.mat', 'genes', 'sconnDen', 'sconnLen', ...
        'ROIsize', 'emp_atrophy');
    
    N_regions = 41;
    v = 1;
    dt = 0.01;
    T_total = 10000;
    init_number = 1;
    syn_control = ROIsize;
    prob_stay = 0.5;
    trans_rate = 1;
    % init seed to anterior hip
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

    % Iteratively call SIR simulator and SIR atrophy for each gene pair.
    %
    % SIR simulator returns (1) normal protein expression over regions over
    % time, and (2) misfolded protein expression over regions over time.
    %
    % From this we derive theortical atrophy per region for each gene and
    % correlate with empirical atrophy per region per gene.

    gene_corrs = table('Size', [width(risk_genes)*width(clear_genes), 6], ...
    'VariableTypes', {'string', 'string', 'double', 'double', 'cell', 'cell'}, ...
    'VariableNames', {'risk gene', 'clearance gene', 'correlation', 't', 'simulated atrophy', 'error'});

    index = 1;
    for risk_gene = 1:width(risk_genes)
        % tic
        risk_name = risk_genes.Properties.VariableNames{risk_gene};
        for clear_gene = 1:width(clear_genes)
            clear_name = clear_genes.Properties.VariableNames{clear_gene};
            if strcmp(clear_name, risk_name)
                continue
            end
            %%%%% simulation ------ >>>
            [Rnor_all, Rmis_all] = SIRsimulator(N_regions, v, dt, T_total, ...
                clear_genes(:,clear_gene), risk_genes(:,risk_gene), ...
                sconnLen, sconnDen, ROIsize, seed, syn_control, ...
                init_number, prob_stay, trans_rate);
            [sim_atrophy] = SIRatrophy(Rnor_all, Rmis_all, sconnDen, ...
                N_regions, dt, emp_atrophy, risk_name, clear_name, plotting);
            gene_corrs(index,:) = {risk_name, clear_name, max_corr, tstep, sim_atrophy(:, tstep), err};
            index = index + 1;
        end
        % toc
    end
    
    gene_corrs = sortrows(gene_corrs, 3);
    
    tail(gene_corrs)
    % writetable(gene_corrs, 'results/gene_corrs.csv')
    % writetable(gene_corrs, 'results/gene_corrs.csv', 'WriteMode', 'Append')

end