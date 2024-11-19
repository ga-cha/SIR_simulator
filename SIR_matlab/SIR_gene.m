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

        % intermediate products
        Rnor_all;
        Rmis_all;
        sim_atrophy;

        % output table
        gene_corr;
    end

    methods
        function self = SIR_gene(genes, idx)
            % Rearrange index to  risk and clear gene indices
            [i, j] = ind2sub([genes.n_risk, genes.n_clear], idx);
            % Sets name and expression array of this specific gene pair
            self.risk_name = genes.risk_names.Value(i);
            self.risk_gene = genes.risk_genes.Value(:, i);
            self.clear_name = genes.clear_names.Value(j);
            self.clear_gene = genes.clear_genes.Value(:, j);

            self = init_corr_table(self);
        end

        function self = run_sim(self, params)
            if strcmp(self.risk_name, self.clear_name); return; end

            % sir_simulator and sir_atrophy have extra visualisations,
            % infrequently used and called separately to vis
            [self.Rnor_all, self.Rmis_all] = sir_simulator(params, self, false);
            self = sir_atrophy(params, self, false);
            [bgs, cobre, hcpep, stages, tstep] = sir_corr(params, self);
            avg_max = mean([bgs, cobre, hcpep]);

            self.gene_corr = {self.risk_name, self.clear_name,              ...
                avg_max, bgs, cobre, hcpep, stages, tstep};
        end

        function self = init_corr_table(self)
            self.gene_corr = table(                                         ...
                'Size', [1, 8],                                             ...
                'VariableTypes', {'string', 'string', 'double', 'double',   ...
                'double', 'double', 'double', 'double'},                    ...
                'VariableNames', {'risk gene', 'clearance gene',            ...
                'correlation', 'bgs correlation', 'cobre correlation',      ...
                'hcpep correlation', 'stages correlation', 'bgs t'});
        end
    end
end