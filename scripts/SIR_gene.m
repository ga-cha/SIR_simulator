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
        idx;

        % output correlations
        max_corr;               % ([site_name, site_max, site_t] x sites)
        gene_corr = {};
        col_names = {'risk gene', 'clearance gene', 'seed', 't_max', 't_05'};
        t_max = 0;
        t_05 = 0;

    end

    methods
        function self = SIR_gene(params, idx)
            if nargin == 0; return; end % needs empty constructor for parfor

            self.idx = idx;
            % Rearrange index to  risk and clear gene indices
            [i, j] = ind2sub([params.n_risk, params.n_clear], idx);
            % Sets name and expression array of this specific gene pair
            self.risk_name = params.risk_names.Value(i);
            self.risk_gene = params.risk_genes.Value(:, i);
            self.clear_name = params.clear_names.Value(j);
            self.clear_gene = params.clear_genes.Value(:, j);

            self.max_corr = table('Size', [params.n_sites, 3], ...
                'VariableTypes', {'string', 'double', 'double'}, ...
                'VariableNames', {'site_name', 'site_corr', 'site_t'});

        end

        function self = run_gene(self, params)
            if strcmp(self.risk_name, self.clear_name); return; end

            exp_gene = run_model(self, params, 1);

            if params.null == "none"
                self = exp_gene;
            else
                params = params.set_null();
                if params.null == "spatial"
                    self = run_spatial(self, params);
                elseif params.null == "rewired"
                    self = run_rewired(self, params);
                elseif params.null == "gene"
                    self = run_gene_null(self, params);
                end
                % Calculate p value of experimental correlation from null distribution
                % correlation column is the first column after preamble
                null_corrs = cell2mat(self.gene_corr(:, width(self.col_names)+1));
                exp_corr = exp_gene.gene_corr{:, width(exp_gene.col_names)+1};
                p = tail_approx(self, null_corrs, exp_corr, params.vis);

                % Assign null summary statistics to output table
                % self.gene_corr = [self.gene_corr(1, 1:2), params.seed, p];
                self.gene_corr = [exp_gene.gene_corr, p];

                if params.vis; plot_null(self, null_corrs, exp_corr, p); end
            end

            % Cleanup
            self.max_corr = table();
        end

        function p = tail_approx(self, null_corrs, exp_corr, vis)
            % Use palm_pareto (Winkler et al.) for tail/GPD approximation.
            % https://github.com/andersonwinkler/PALM
            null_corrs = double(null_corrs(:));
            assert(~isempty(null_corrs), 'Null correlations are empty');

            try
                [p, apar, kpar, upar] = palm_pareto(exp_corr, null_corrs, false, 0.1, false);
            catch
                p = palm_pareto(exp_corr, null_corrs, false, 0.1, false);
            end

            % Visualization: reconstruct tail/excess for plotting if possible
            if vis && exist('upar','var') && ~isempty(upar) && exp_corr > upar
                tail_data = null_corrs(null_corrs > upar);
                if ~isempty(tail_data)
                    excesses = tail_data - upar;
                    excess_exp = exp_corr - upar;
                    parmhat = [kpar, apar]; % [shape, scale] matching plot_gpd expectation
                    plot_gpd(self, excesses, excess_exp, parmhat);
                end
            end
        end
        
        function self = calculate_correlations(self, params, sim_atrophy, i)
            self = sir_corr(params, self, sim_atrophy);
            max_corrs = self.max_corr{:, 2}';
            % max_t = self.max_corr{1, 3};
            self.gene_corr(i, :) =                                                 ...
                [{self.risk_name, self.clear_name, params.seed, self.t_max, self.t_05},     ...
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
                params.emp_atr = params.null_atr{i};
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

        function self = run_gene_null(self, params)
            % run_gene_null produces a different simulated result to run_model
            % by setting different gene expression profiles
            % this simulation is correlated against empirical atrophy
            for i = 1:1000
                params = set_gene_null(params, i);

                % TODO: refactor
                [j, k] = ind2sub([params.n_risk, params.n_clear], self.idx);
                self.risk_gene = params.risk_genes.Value(:, j);
                self.clear_gene = params.clear_genes.Value(:, k);

                self = self.run_model(params, i);
            end
        end

        function [] = plot_null(self, null_corrs, exp_corr, p)
            % Set box plot whiskers to 95% CI. Taken from:
            % https://www.mathworks.com/matlabcentral/answers/171414-how-to-show-95-quanile-in-a-boxplot
            q95=norminv(0.95);
            q3=norminv(.75);
            w95=(q95-q3)/(2*q3);

            figure('Color','white', 'Position', [100 100 350 500]);
            hold on;
            boxplot(null_corrs, 'whisker', w95)

            s=swarmchart(ones(size(null_corrs)), null_corrs, 5, 'MarkerEdgeAlpha', 0.2);
            s.XJitterWidth = 0.15;
            plot(1, exp_corr, '_', 'MarkerSize', 35, 'LineWidth', 2);

            t = title({['Pearson correlation = ', num2str(exp_corr)], ['p = ', num2str(p)]});
            t.FontWeight = 'normal';
            xlabel({"risk gene: " + self.risk_name, "clearance gene: " + self.clear_name})
            ylabel('correlation')
            xlim([0.8 1.2])
            ylim([min([null_corrs; exp_corr])-0.1, max([null_corrs; exp_corr])+0.1])
            % ylim([-0.2 0.8])

            hold off;
        end

        function [] = plot_gpd(~, excesses, excess_exp, parmhat)
            % Plot GPD fit to null_corrs > threshold and exp_corr
            figure('Color','white');
            hold on;
            % Histogram of tail data
            histogram(excesses, 'Normalization', 'pdf', 'FaceAlpha', 0.5);
            % GPD fit curve
            x_fit = linspace(0, max([excesses; excess_exp]), 100);
            y_fit = gppdf(x_fit, parmhat(1), parmhat(2));
            plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
            % plot exp_corr
            xline(excess_exp, '--', 'LineWidth', 2);
            xlabel('Excess above threshold');
            ylabel('PDF');
            legend({'Tail data', 'GPD fit', 'exp\_corr'}, 'Location', 'best');
            title('GPD fit to null tail');
            hold off;
        end
    end
end