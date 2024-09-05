% gene tables are indexed as T(:, ["gene1", "gene2"])
clear = ["RFWD3"];
opt.risk = ["PARVA"];
opt.parc = "S132";
opt.dbg = false;
opt.vis = false;
opt.null = "none";

ws = 'data/workspace_' + opt.parc + '.mat';
% load workspace variables in parameter object 
params = SIRparameters(opt);
params = set_atrophy(params, ws);
params = set_netw(params, ws);
ws_null = 'data/workspace_' + opt.parc + '_null.mat';
params = set_null(params, opt, ws_null);
% load gene expression data in genes parameter object
genes = SIRgenes(opt, clear, ws);
gene = SIRgene(genes, 1);

%% Experimental gene correlation
[genes] = SIRiterator(params, genes);
genes.gene_corrs = rmmissing(genes.gene_corrs);
exp_corr = genes.gene_corrs.correlation;

%% Null gene correlation
params.null = "spatial";
opt.null = "spatial";
params = set_null(params, opt, ws_null);

genes = SIRgenes(opt, clear, ws);
[genes] = SIRiterator(params, genes);

genes.gene_corrs = rmmissing(genes.gene_corrs);
null_corrs = genes.gene_corrs.correlation;

%% t-test and plot
[h,p]=ttest2(exp_corr,null_corrs);

% Set box plot whiskers to 95% CI. Taken from:
% https://www.mathworks.com/matlabcentral/answers/171414-how-to-show-95-quanile-in-a-boxplot
q95=norminv(0.95);
q3=norminv(.75);
w95=(q95-q3)/(2*q3);

figure;
hold on
boxplot(null_corrs, 'whisker', 0.7193)
% s=swarmchart(ones(length(null_corrs),1),null_corrs,5,'MarkerEdgeAlpha',0.5);
s=swarmchart(ones(length(null_corrs),1),null_corrs,5,'MarkerEdgeAlpha',0.2);
s.XJitterWidth=0.15;
ylim([0 0.85])
plot(1,exp_corr,'_','LineWidth', 1.2,'MarkerSize',35)
hold off

t = title({['Spearman correlation = ', num2str(exp_corr)], ['p = ', num2str(p)]});
t.FontWeight = 'normal';
xlabel({"risk gene: " + gene.risk_name, "clearance gene: " + gene.clear_name})
ylabel('correlation')
xlim([0.8 1.2])
ylim([-0.2 0.7])

function [genes] = SIRiterator(params, genes)
    n = genes.n_risk * genes.n_clear;
    for idx = 1:n; genes = iterate(idx, params, genes); end
    
    assert (~isempty(genes.gene_corrs));
    genes.gene_corrs = rmmissing(genes.gene_corrs);
end

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