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

    ws = 'data/workspace_' + opt.parc + '.mat';
    % load workspace variables in parameter object 
    params = SIRparameters(opt);
    params = set_atrophy(params, ws);
    params = set_netw(params, ws);
    ws_null = 'data/workspace_' + opt.parc + '_null.mat';
    params = set_null(params, opt, ws_null);
    % load gene expression data in genes parameter object
    genes = SIRgenes(opt, clear, ws);

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
tic
    gene = SIRgene(genes, idx);
    if strcmp(gene.clear_name, gene.risk_name); return; end

    % SIRsimulator and SIRatrophy have visualization as a separate
    % parameter, switched off.
    [gene.Rnor_all, gene.Rmis_all] = SIRsimulator(params, gene, false);
    gene = SIRatrophy(params, gene, false);
toc
tic
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
toc
end



% % First we generate simulated atrophy 
% [Rnor_all, Rmis_all] = SIRsimulator(N_regions, v, dt, T_total, ...
%     clear_gene, risk_gene, sconnLen, sconnDen, ROIsize, seed, ...
%     syn_control, init_number, prob_stay, trans_rate);
% [sim_atrophy] = SIRatrophy(Rnor_all, Rmis_all, sconnDen, N_regions, dt);
% 
% % Then we determine the correlation with null atrophy. Null atrophy is an 
% % input of scrambled empirical atrophy, preserving spatial autocorrelation
% 
% nulls = 1000;
% null_corrs = zeros(nulls,1);
% 
% for i = 1:nulls
%     [null_corrs(i), ~] = SIRcorr(sim_atrophy, null_atrophy(i, :), T_total);
% end
% 
% % experimental atrophy correlation
% [exp_corr, ~] = SIRcorr(sim_atrophy, emp_atrophy, T_total);
% 
% % p-value
% [h,p]=ttest2(exp_corr,null_corrs);
% disp(['p = ', num2str(p)])
% 
% % Set box plot whiskers to 95% CI. Taken from:
% % https://www.mathworks.com/matlabcentral/answers/171414-how-to-show-95-quanile-in-a-boxplot
% q95=norminv(0.95);
% q3=norminv(.75);
% w95=(q95-q3)/(2*q3);
% 
% figure;
% hold on
% boxplot(null_corrs, 'whisker', 0.7193)
% swarmchart(ones(length(null_corrs),1),null_corrs,5,'MarkerEdgeAlpha',0.5,'jitter','on')
% ylim([0 0.85])
% plot(1,exp_corr,'_','LineWidth', 1.2,'MarkerSize',35)
% hold off