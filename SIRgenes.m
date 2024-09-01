% SIRgenes.m
%
% Gabriella Chan 25/06/24
% gabriella.chan@monash.edu
%
% Takes a gene expression table, a string array of clearance gene names,
% and an optional string array of risk gene names.
%
%

classdef SIRgenes
    properties
        risk_names;
        clear_names;
        clear_genes;             % clear_gene: clear_gene gene expression (zscore, N_regions * 1 vector) (empirical clear_gene expression)
        risk_genes;              % risk_gene: risk_gene gene expression after normalization (zscore, N_regions * 1 vector) (empirical risk_gene expression)
        n_risk;
        n_clear;
        gene_corrs;
    end

    methods
        % Sets risk and clearance gene arrays
        function self = SIRgenes(opt, clear_names, ws)
            load(ws, 'gene_expr');
            % TODO: get this function working
            risk_names = set_risk_names([gene_expr.Properties.VariableNames], opt);
            % if ~isfield(opt, 'risk')
            %     if opt.parc == "S132_PCA"               % 65 risk PCs
            %         opt.risk = gene_expr.Properties.VariableNames(1:65);
            %     elseif opt.parc == "S332_PCA"           % 165 risk PCs
            %         opt.risk = gene_expr.Properties.VariableNames(1:165);
            %     else
            %         opt.risk = string(gene_expr.Properties.VariableNames);
            %     end
            % end
            % risk_names = opt.risk;

            % mask input genes that exist in gene expression table
            % then generate gene expression table of masked genes only
            clear_mask = ismember(clear_names, gene_expr.Properties.VariableNames);
            self.clear_names = clear_names(clear_mask);
            self.clear_genes = gene_expr(:, clear_names);
            self.n_clear = width(self.clear_genes);

            risk_mask = ismember(risk_names, gene_expr.Properties.VariableNames);
            self.risk_names = risk_names(risk_mask);
            self.risk_genes = gene_expr(:, risk_names);
            self.n_risk = width(self.risk_genes);

            % preslice variables for parfor loop
            self = preslice(self);

            % initialise gene correlation output table
            self = init_gene_corrs(self, opt);

            assert (~isempty(self.clear_genes) && ~isempty(self.risk_genes), ...
                "No valid gene combinations");
        end

        function risk_names = set_risk_names(gene_names, opt)
            if ~isfield(opt, 'risk')
                if opt.parc == "S132_PCA"               % 65 risk PCs
                    opt.risk = gene_names(1:65);
                elseif opt.parc == "S332_PCA"           % 165 risk PCs
                    opt.risk = gene_names(1:165);
                else
                    opt.risk = string(gene_names);
                end
            end
            risk_names = opt.risk;
        end

        function self = init_gene_corrs(self, opt)
            dim = 1;
            if opt.null == "spatial"; dim = 1000; end
            self.gene_corrs = table(                                                ...
                'Size', [width(self.risk_genes) * width(self.clear_genes) * dim, 8],      ...
                'VariableTypes', {'string', 'string', 'double', 'double', 'double', ...
                'double', 'double', 'double'},                                      ...
                'VariableNames', {'risk gene', 'clearance gene', 'correlation',     ...
                'bgs correlation', 'cobre correlation', 'hcpep correlation',        ...
                'stages correlation', 'bgs t'});
        end

        function self = preslice(self)
            self.risk_names = parallel.pool.Constant(self.risk_names);
            self.clear_names = parallel.pool.Constant(self.clear_names);
            self.risk_genes = parallel.pool.Constant(table2array(self.risk_genes));
            self.clear_genes = parallel.pool.Constant(table2array(self.clear_genes));
        end
    end
end