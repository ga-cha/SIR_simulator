% gene tables are indexed as T(:, ["gene1", "gene2"])
% clear_names = ["RNF144A"];
% opt.risk_names = ["CDK19"];
clear_names = ["gene350"];
opt.risk_names = ["gene2453"];
% opt.parc = "S332";
opt.parc = "S332_rand";
opt.dbg = false;
opt.vis = false;
opt.null = "none";

% load workspace variables in parameter object 
params = SIR_parameters(opt);

% load gene expression data in genes parameter object
genes = SIR_genes(opt, clear_names);

% Experimental gene correlation
gene_corrs = run_sim(params, genes);
exp_corr = gene_corrs.correlation;

% Null gene correlation
params.null = "spatial";
opt.null = "spatial";
params = params.set_null(opt);

gene_corrs = run_sim(params, genes);
null_corrs = gene_corrs.correlation;

% t-test and plot
[h,p]=ttest2(exp_corr,null_corrs);

% Set box plot whiskers to 95% CI. Taken from:
% https://www.mathworks.com/matlabcentral/answers/171414-how-to-show-95-quanile-in-a-boxplot
q95=norminv(0.95);
q3=norminv(.75);
w95=(q95-q3)/(2*q3);

figure;
hold on;
boxplot(null_corrs, 'whisker', w95)
s=swarmchart(ones(length(null_corrs),1),null_corrs,5,'MarkerEdgeAlpha',0.2);
s.XJitterWidth=0.15;
ylim([0 0.85])
plot(1,exp_corr,'_','LineWidth', 1.2,'MarkerSize',35)

t = title({['Pearson correlation = ', num2str(exp_corr)], ['p = ', num2str(p)]});
t.FontWeight = 'normal';
xlabel({"risk gene: " + opt.risk_names, "clearance gene: " + clear_names})
ylabel('correlation')
xlim([0.8 1.2])
ylim([-0.2 0.7])

hold off;

set(gcf, 'Position',[100 100 350 500])


function gene_corrs = run_sim(params, genes)
    n = genes.n_risk * genes.n_clear;
    gene_objs = SIR_gene.empty(n, 0);

    for idx = 1:n
        gene = SIR_gene(genes, idx);
        if params.null == "spatial"
            gene_objs(idx) = gene.null_sim(params);
        else 
            gene_objs(idx) = gene.run_sim(params);
        end
    end
    
    gene_corrs = assemble_gene_corrs(gene_objs);
end

% assemble table from gene pairs, preallocating for speed
function [gene_corrs] = assemble_gene_corrs(gene_objs)
    % collate correlation from each gene object into a single cell array
    gene_corr_cells = cell(length(gene_objs), 1);
    parfor i = 1:length(gene_objs)
        if isprop(gene_objs(i), 'gene_corr')
            gene_corr_cells{i} = gene_objs(i).gene_corr;
        end
    end
    % concatenate cell array into a single table
    C = vertcat(gene_corr_cells{:});
    assert (~isempty(C));
    gene_corrs = cell2table(C, 'VariableNames', {'risk gene',               ...
        'clearance gene', 'correlation', 'bgs correlation',                 ...
        'cobre correlation', 'hcpep correlation', 'stages correlation'});
end