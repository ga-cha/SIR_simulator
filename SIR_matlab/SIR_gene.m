% SIRgene.m
%
% Gabriella Chan 17/08/24
% gabriella.chan@monash.edu
%
% Holds a risk/clear gene pair, intermediate values, and outputs
%

classdef SIR_gene
    properties
        risk_name;
        risk_gene;
        clear_name;
        clear_gene;

        % output correlations
        max_atr = table();
        max_corr = table();
        gene_corr = {};
    end

    methods
        function self = SIR_gene(genes, idx)
            if nargin == 0; return; end % needs empty constructor for parfor

            % Rearrange index to  risk and clear gene indices
            [i, j] = ind2sub([genes.n_risk, genes.n_clear], idx);
            % Sets name and expression array of this specific gene pair
            self.risk_name = genes.risk_names.Value(i);
            self.risk_gene = genes.risk_genes.Value(:, i);
            self.clear_name = genes.clear_names.Value(j);
            self.clear_gene = genes.clear_genes.Value(:, j);
        end

        % Simulate normal and misfolded protein growth
        % Calculate simuated atrophy over time
        % Calculate correlations over time
        function self = run_sim(self, params)
            if strcmp(self.risk_name, self.clear_name); return; end
            
            % sir_simulator and sir_atrophy have extra visualisations,
            % infrequently used and called separately to vis
            proteins = sir_simulator(params, self, false);
            sim_atrophy = sir_atrophy(params, proteins, false);

            self = sir_corr(params, self, sim_atrophy);

            % mean_max_corr = mean(table2array(self.max_corr(:, 2)));
            max_corrs = table2array(self.max_corr(:, 2))';
    
            self.gene_corr = [self.gene_corr;                           ...
                {self.risk_name, self.clear_name},       ...
                num2cell(max_corrs)];
        end

        % spatial null replicates a lot of above
        % TODO: refactor
        function self = null_sim(self, params)
            if strcmp(self.risk_name, self.clear_name); return; end
            
            proteins = sir_simulator(params, self, false);
            sim_atrophy = sir_atrophy(params, proteins, false);

            for i = 1:1000
                self.max_atr = table();
                self.max_corr = table();

                self = sir_corr(params, self, sim_atrophy, i);
    
                % mean_max_corr = mean(table2array(self.max_corr(:, 2)));
                max_corrs = table2array(self.max_corr(:, 2))';
        
                self.gene_corr = [self.gene_corr;                           ...
                    {self.risk_name, self.clear_name},       ...
                    num2cell(max_corrs)];
            end
        end

    end
end