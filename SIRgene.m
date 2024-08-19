% SIRgene.m
%
% Gabriella Chan 17/08/24
% gabriella.chan@monash.edu
%
% the gene struct holds a risk/clear gene pair and intermediate outputs
%

classdef SIRgene
    properties
        risk_name;
        risk_gene;
        clear_name;
        clear_gene;
    end

    methods
        % Sets risk and clearance gene arrays
        function self = SIRgene(genes, idx)
            % Rearrange index to  risk and clear gene indices
            [i, j] = ind2sub([genes.n_risk, genes.n_clear], idx);

            self.risk_name = genes.risk_names.Value(i);
            self.risk_gene = genes.risk_genes.Value(:, i);
            self.clear_name = genes.clear_names.Value(j);
            self.clear_gene = genes.clear_genes.Value(:, j);
        end
    end
end