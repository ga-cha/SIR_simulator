% SIRgenes.m
%
% Gabriella Chan 25/06/24
% gabriella.chan@monash.edu
%
% Takes a gene expression table, a string array of clearance gene names,
% and an optional string array of risk gene names.
%
%

classdef SIR_genes
    properties
        risk_names;
        clear_names;
        clear_genes;             % clear_gene: clear_gene gene expression (zscore, N_regions * 1 vector) (empirical clear_gene expression)
        risk_genes;              % risk_gene: risk_gene gene expression after normalization (zscore, N_regions * 1 vector) (empirical risk_gene expression)
        n_risk;
        n_clear;
    end

    methods
        % Sets risk and clearance gene arrays
        function self = SIR_genes(opt, clear_names)
            ws = 'data/workspace_' + opt.parc + '.mat';
            load(ws, 'gene_expr');

            risk_names = set_risk_names(self, opt, gene_expr);

            % mask out gene names that don't exist
            % then generate gene expression table of masked genes only
            clear_mask = ismember(clear_names, gene_expr.Properties.VariableNames);
            self.clear_names = clear_names(clear_mask);
            self.clear_genes = gene_expr(:, self.clear_names);
            self.n_clear = width(self.clear_names);

            risk_mask = ismember(risk_names, gene_expr.Properties.VariableNames);
            self.risk_names = risk_names(risk_mask);
            self.risk_genes = gene_expr(:, self.risk_names);
            self.n_risk = width(self.risk_names);

            % preslice variables for parfor loop
            self = preslice(self);

            % TODO: check this break statement
            assert (~isempty(self.clear_genes) && ~isempty(self.risk_genes), ...
                "No valid gene combinations");
        end

        function risk_names = set_risk_names(self, opt, gene_expr)          %#ok<INUSD>
            % self looks unnecessary, but is for some reason necessary

            if ~isfield(opt, 'risk_names')
                if opt.parc == "S132_PCA"               % 65 risk PCs
                    opt.risk_names = gene_expr.Properties.VariableNames(1:65);
                elseif opt.parc == "S332_PCA"           % 165 risk PCs
                    opt.risk_names = gene_expr.Properties.VariableNames(1:165);
                else
                    opt.risk_names = string(gene_expr.Properties.VariableNames);
                end
            end

            % opt.risk_names = ["LAMP5", "PARVA"];
            % opt.risk_names = ["gene5", "gene6"];
            risk_names = opt.risk_names;

            assert(opt.vis == false | (length(risk_names) < 5 & opt.null == "none"), ...
                "It is not recommended to visualise this many combinations")
        end

        function self = preslice(self)
            self.risk_names = parallel.pool.Constant(self.risk_names);
            self.clear_names = parallel.pool.Constant(self.clear_names);
            self.risk_genes = parallel.pool.Constant(table2array(self.risk_genes));
            self.clear_genes = parallel.pool.Constant(table2array(self.clear_genes));
        end        
    end
end