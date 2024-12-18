% main_null_rewired.m
%                                                                                                                                    
% Package based off https://github.com/yingqiuz/SIR_simulator
% Zheng, Ying-Qiu, et al. PLoS biol. 17.11 (2019): e3000495.
%
% Gabriella Chan 13/08/23
% gabriella.chan@monash.edu
% Monash University
%
% We generate simulated atrophy of our gene pair with each different
% rewired null, and record the correlation with empirical atrophy. We then
% compute the correlation of simulated atrophy with empirical atrophy when
% using the true structural connectome density.
% We calculate the p-value of experimental correlation against nulls, and 
% plot on a box-and-whiskers plot.


function [] = main_null_rewired(clear_genes_in, risk_genes_in, debuggy, plotting)


    % load gene expressions, real atrophy, ROIsize, functional connectivity...
    % simulated structural connectivity lengths are described in greater 
    % detail in Zheng PLoS biol. (2019)
    % The implementation is given in SIR_utils/null_ROI_dist.m
    load('data/workspace_DK_pos_atrophy.mat', 'genes', 'sconnDen', 'sconnLen', ...
        'ROIsize', 'emp_atrophy_bgs', 'emp_atrophy_cobre', 'emp_atrophy_hcpep', ...
        'emp_atrophy_stages', 'risk_names');
    % Load simulated structural connectivity density. Derived in
    % null_rewire.m, from BCT implementation of Maslov Sneppen
    load('data/rewire/sconnDen_sim.mat', 'sconnDen_sim');
    
    N_regions = 41;
    v = 1;
    dt = 0.01;
    T_total = 10000;
    init_number = 1;
    syn_control = ROIsize;
    prob_stay = 0.5;
    trans_rate = 1;
    plotting = false;
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
    % single risk/clearance gene pair input, as tables
    clear_gene = genes(:, 'REEP4');
    risk_gene = genes(:, 'PARVA');
    nulls = 1000;
    null_corrs = zeros(nulls,1);
    
    for i = 1:nulls
        [gene_corrs] = SIRiterator(N_regions, v, dt, T_total, ...
            clear_gene, risk_gene, sconnLen, sconnDen_sim(:,:,i), ROIsize, seed, ...
            init_number, prob_stay, trans_rate, emp_atrophy_bgs, ...
        emp_atrophy_cobre, emp_atrophy_hcpep, emp_atrophy_stages, plotting);
        null_corrs(i) = gene_corrs.correlation;
    end
    
    % experimental atrophy
    [gene_corrs] = SIRiterator(N_regions, v, dt, T_total, ...
        clear_gene, risk_gene, sconnLen, sconnDen, ROIsize, seed, ...
        init_number, prob_stay, trans_rate, emp_atrophy_bgs, ...
        emp_atrophy_cobre, emp_atrophy_hcpep, emp_atrophy_stages, plotting);
    exp_corr = gene_corrs.correlation;
    % p-value
    [h,p]=ttest2(exp_corr,null_corrs);
    disp(['p = ', num2str(p)])
    
    % Set box plot whiskers to 95% CI. Taken from:
    % https://www.mathworks.com/matlabcentral/answers/171414-how-to-show-95-quanile-in-a-boxplot
    q95=norminv(0.95);
    q3=norminv(.75);
    w95=(q95-q3)/(2*q3);
    
    figure;
    hold on
    boxplot(null_corrs, 'whisker', 0.7193)
    swarmchart(ones(length(null_corrs),1),null_corrs,5,'MarkerEdgeAlpha',0.5,'jitter','on')
    ylim([0 0.85])
    plot(1,exp_corr,'_','LineWidth', 1.2,'MarkerSize',35)
    hold off

end