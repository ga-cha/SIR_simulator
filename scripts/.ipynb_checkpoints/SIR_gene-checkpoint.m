% SIRgene.m
%
% Gabriella Chan 17/08/24
% gabriella.chan@monash.edu
%
% Holds a risk/clear gene pair and outputs
%

classdef SIR_gene
    properties
        risk_name;
        risk_gene;
        clear_name;
        clear_gene;

        % output correlations
        max_corr;               % ([site_name, site_max, site_t] x sites)
        gene_corr = {};
        col_names = {'risk gene', 'clearance gene', 'seed', 't_90', 't_05'};
        t_90 = 0;
        t_05 = 0;

    end

    methods
        function self = SIR_gene(genes, params, idx)
            if nargin == 0; return; end % needs empty constructor for parfor

            % Rearrange index to  risk and clear gene indices
            [i, j] = ind2sub([genes.n_risk, genes.n_clear], idx);
            % Sets name and expression array of this specific gene pair
            self.risk_name = genes.risk_names.Value(i);
            self.risk_gene = genes.risk_genes.Value(:, i);
            self.clear_name = genes.clear_names.Value(j);
            self.clear_gene = genes.clear_genes.Value(:, j);

            self.max_corr = table('Size', [params.n_sites, 3], ...
                'VariableTypes', {'string', 'double', 'double'}, ...
                'VariableNames', {'site_name', 'site_corr', 'site_t'});

        end

        function self = run_gene(self, params)
            if strcmp(self.risk_name, self.clear_name); return; end

            if params.null == "none"
                self = run_model(self, params, 1);
            else
                if params.null == "spatial"
                    self = run_spatial(self, params);
                elseif params.null == "rewired"
                    self = run_rewired(self, params);
                end
                % Assign null summary statistics to output table
                % correlation column is the first column after preamble
                corrs = cell2mat(self.gene_corr(:, width(self.col_names)+1));
                self.gene_corr = [self.gene_corr(1, 1:2), params.seed,      ...
                    mean(corrs), std(corrs)];
            end

            % Cleanup
            self.max_corr = table();
        end

        function self = calculate_correlations(self, params, sim_atrophy, i)
            self = sir_corr(params, self, sim_atrophy);
            max_corrs = self.max_corr{:, 2}';
            % max_t = self.max_corr{1, 3};
            self.gene_corr(i, :) =                                                 ...
                [{self.risk_name, self.clear_name, params.seed, self.t_90, self.t_05},     ...
                num2cell(max_corrs)];
        end

        function self = run_model(self, params, i)
            % Usage: i is row index of output gene correlation table

            % sir_simulator and sir_atrophy have extra visualisations,
            % infrequently used and called separately to vis
            
            % Simulate normal and misfolded protein growth
            proteins = sir_simulator(params, self, false);

            % Calculate simuated atrophy over time
            sim_atrophy = sir_atrophy(params, proteins, false);
            % Calculate correlations over time
            self = calculate_correlations(self, params, sim_atrophy, i);
        end
        
        function self = run_spatial(self, params)
            % run_spatial produces the same simulated result as run_model
            % this simulation is correlated against null maps of atrophy
            proteins = sir_simulator(params, self, false);
            sim_atrophy = sir_atrophy(params, proteins, false);
            for i = 1:1000
                params = set_spatial(params, i); 
                self = calculate_correlations(self, params, sim_atrophy, i);
            end
        end 

        function self = run_rewired(self, params)
            % run_rewired produces a different simulated result to run_model
            % by setting different connectome weights
            % this simulation is correlated against empirical atrophy
            for i = 1:1000
                params.sc_weight = params.null_weight(:, :, i);
                self = self.run_model(params, i);
            end
        end
    end
end