% SIRiterator
%
% Gabriella Chan 28/06/23
% gabriella.chan@monash.edu
% Monash University
%
%
% This module serves as an intermediate layer to call SIR simulator
% and SIR atrophy.
%
% We iterate through all gene pairs and their gene expression and call SIR
% simulator to store an array of normal protein expression over regions 
% over time, and misfolded protein expression over regions over time.
%
% From this we derive theortical atrophy per region for each gene and
% correlate with empirical atrophy per region per gene.

function [gene_corrs, index] = SIRiterator(gene_corrs, index, col3, N_regions, v, dt, ...
    T_total, clear_genes, risk_genes, sconnLen, sconnDen, ROIsize, ...
    seed, init_number, prob_stay, trans_rate, emp_atrophy_bgs, ...
    emp_atrophy_cobre, emp_atrophy_hcpep, emp_atrophy_stages, plotting)

    for risk_gene = 1:width(risk_genes)
        % tic
        risk_name = risk_genes.Properties.VariableNames{risk_gene};
        for clear_gene = 1:width(clear_genes)
            clear_name = clear_genes.Properties.VariableNames{clear_gene};
            if strcmp(clear_name, risk_name)
                continue
            end
            %%%%% simulation ------ >>>
            [Rnor_all, Rmis_all] = SIRsimulator(N_regions, v, dt, T_total, ...
                clear_genes(:,clear_gene), risk_genes(:,risk_gene), sconnLen, ...
                sconnDen, ROIsize, seed, init_number, prob_stay, trans_rate, plotting);
            [sim_atrophy] = SIRatrophy(Rnor_all, Rmis_all, ROIsize, sconnDen, ...
                N_regions, dt, plotting);
            [bgs_max, cobre_max, hcpep_max, stages_max, tstep, ~] = SIRcorr(sim_atrophy, emp_atrophy_bgs, ...
                emp_atrophy_cobre, emp_atrophy_hcpep, emp_atrophy_stages, ...
                risk_name, clear_name, plotting);
            % avg_max = mean(cat (1,bgs_max, cobre_max, hcpep_max, stages_max));
            end_max = mean(cat (1,bgs_max, cobre_max));
            % Note col3 is the test variable; e.g. seed, null #, etc
            gene_corrs(index,:) = {risk_name, clear_name, col3, end_max, bgs_max, cobre_max, hcpep_max, stages_max, tstep};
            index = index + 1;
        end
        % toc
    end
end
