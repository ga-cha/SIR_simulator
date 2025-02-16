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

        function self = calculate_correlations(self, params, sim_atrophy)
            self = sir_corr(params, self, sim_atrophy);
            max_corrs = table2array(self.max_corr(:, 2))';
            self.gene_corr = [self.gene_corr;                           ...
                {self.risk_name, self.clear_name}, num2cell(max_corrs)];
        end

        function self = run_gene(self, params)
            if strcmp(self.risk_name, self.clear_name); return; end

            if params.null == "none"
                self = run_model(self, params);
            elseif params.null == "spatial"
                self = run_spatial(self, params);
            elseif params.null == "rewired"
                self = run_rewired(self, params);
            end

            % Cleanup
            self.max_atr = table();
            self.max_corr = table();
        end

        function self = run_model(self, params)
            % sir_simulator and sir_atrophy have extra visualisations,
            % infrequently used and called separately to vis
            
            % Simulate normal and misfolded protein growth
            proteins = sir_simulator(params, self, false);
            % Calculate simuated atrophy over time
            sim_atrophy = sir_atrophy(params, proteins, false);
            % Calculate correlations over time
            self = calculate_correlations(self, params, sim_atrophy);
        end
        
        function self = run_spatial(self, params)
            % run_spatial produces the same simulated result as run_model
            % this simulation is correlated against null maps of atrophy
            proteins = sir_simulator(params, self, false);
            sim_atrophy = sir_atrophy(params, proteins, false);
            for i = 1:1000
                params = set_spatial(params, i); 
                self = calculate_correlations(self, params, sim_atrophy);
            end
        end 

        function self = run_rewired(self, params)
            % run_rewired produces a different simulated result to run_model
            % by setting different connectome weights
            % this simulation is correlated against empirical atrophy
            for i = 1:1000
                params.sconnDen = params.null_den(:, :, i);
                self = self.run_model(params);
            end
        end
    end
end