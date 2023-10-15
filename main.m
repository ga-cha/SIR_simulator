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
% a script to simulate atrophy accrual due to the accumulation of misfolded
% protein aggregates

function [] = main(clear_gene)
    % load gene expressions, real atrophy, ROIsize, functional connectivity...
    load('data_gc/GC_workspace.mat', 'genes', 'sconnDen', 'sconnLen', ...
        'ROIsize', 'emp_atrophy');
    
    N_regions = 41;
    v = 1;
    dt = 0.01;
    T_total = 10000;
    init_number = 1;
    syn_control = ROIsize;
    prob_stay = 0.5;
    trans_rate = 1;
    % init seed to hip
    seed = 40;
    
    % clearance and risk genes are input as gene x region expression tables
    try
        clear_genes = genes(:, clear_gene);
    catch 
        return;
    end
    risk_genes = genes;
    
    [gene_corrs, ~] = SIRiterator(N_regions, v, dt, T_total, ...
        clear_genes, risk_genes, sconnLen, sconnDen, ROIsize, seed, ...
        syn_control, init_number, prob_stay, trans_rate, emp_atrophy);
    
    % tail(gene_corrs)
    writetable(gene_corrs, 'results_gc/gene_corrs.csv', 'WriteMode', 'Append')

end