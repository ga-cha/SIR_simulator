% gene tables are indexed as T(:, ["gene1", "gene2"])
clear_names = ["RNF144A"];
opt.risk_names = ["TCAP"];
opt.seed = 38;
opt.parc = "S132";
opt.dbg = false;
opt.vis = false;
opt.pf = true;
opt.init = 1;
opt.null = "none";
opt.beta = 1/sqrt(2);

% load workspace variables in parameter object 
params = SIR_parameters(opt);
p = gcp('nocreate'); 
if (isempty(p) && params.pf == true)
    parpool('Threads'); 
end

% load gene expression data in genes parameter object
genes = SIR_genes(opt, clear_names);

% Experimental gene correlation
gene_corrs = run_sim(params, genes);
exp_corr = gene_corrs.correlation;

% Null gene correlation
opt.null = "spatial";
params = SIR_parameters(opt);

gene_corrs = run_sim(params, genes);
% null_corrs = gene_corrs.correlation;
null_corrs = gene_corrs.diagnosis_beta;

% t-test and plot
% [h,p]=ttest2(exp_corr,null_corrs);

% Fit a kernel density estimator
[f, xi] = ksdensity(null_corrs, 1, 'Function', 'cdf');
p = 1 - f;
if p == 0; p = eps; end

% Fit a kernel density estimator (PDF)
% [f_pdf, xi_pdf] = ksdensity(null_corrs, 'Function', 'pdf');

% Set box plot whiskers to 95% CI. Taken from:
% https://www.mathworks.com/matlabcentral/answers/171414-how-to-show-95-quanile-in-a-boxplot
q95=norminv(0.95);
q3=norminv(.75);
w95=(q95-q3)/(2*q3);

figure('Color','white', 'Position', [100 100 350 500]);
hold on;
boxplot(null_corrs, 'whisker', w95)
% Plot the kernel density estimate
% plot(rescale(f_pdf, 1, 1.15), xi_pdf, 'k-');
s=swarmchart(ones(size(null_corrs)),null_corrs,5,'MarkerEdgeAlpha',0.2);
s.XJitterWidth = 0.15;
plot(1, exp_corr,'_','MarkerSize', 35, 'LineWidth', 2);

t = title({['Pearson correlation = ', num2str(exp_corr)], ['p = ', num2str(p)]});
t.FontWeight = 'normal';
xlabel({"risk gene: " + opt.risk_names, "clearance gene: " + clear_names})
ylabel('correlation')
xlim([0.8 1.2])
% ylim([min([xi_pdf, exp_corr])-0.05, max([xi_pdf, exp_corr])+0.05])
ylim([-0.2 0.8])

hold off;


function gene_corrs = run_sim(params, genes)
    n = genes.n_risk * genes.n_clear;
    gene_objs = SIR_gene.empty(n, 0);

    parfor idx = 1:n
        gene = SIR_gene(genes, params, idx);        
        if params.null == "spatial"
            gene_objs(idx) = gene.run_spatial(params);
        elseif params.null == "rewired"
            gene_objs(idx) = gene.run_rewired(params);
        else
            gene_objs(idx) = gene.run_model(params, 1);
        end
    end

    gene_corrs = assemble_gene_corrs(params, gene_objs);
end

% assemble table from gene pairs, preallocating for speed
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
    % gene_corrs = cell2table(C, "VariableNames", [{'risk gene',               ...
    %         'clearance gene'}, params.emp_atr.Properties.VariableNames]);
    gene_corrs = cell2table(C, "VariableNames", [gene_objs(1).col_names, ...
        params.emp_atr.Properties.VariableNames]);
end