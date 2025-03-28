% SIR_lme.m
%
% Gabriella Chan 23/11/24
% gabriella.chan@monash.edu
%
% Holds results of a linear effects model associated with a single gene
% pair
%

            
            % call from SIR_gene. lme is a class
            % containing an lme model and also some functions used to
            % generate the model and plots, interpretations etc
            % lme = SIR_lme(self.max_atr, params);
            % if params.vis; plot_lme(lme, params); end
            % lme = lme.set_lme();

classdef SIR_lme
    properties
        tbl;
        lme;

        sim_fx;
        lme_corr;
        pval;
    end

    methods
        function self = SIR_lme(sim_atrophy, params)
            % n_rois = params.n_rois;
            n_sites = numel(params.site_names);
                
            % reshape empirical data into a single table
            self.tbl = table();
            for i = 1:n_sites
                site_name = params.site_names{i};
                emp_atr = params.emp_atr(:, site_name);
                emp_atr.Properties.VariableNames = {'EmpAtrophy'};
                sim_atr = sim_atrophy(:, i);
                sim_atr.Properties.VariableNames = {'SimAtrophy'};
                site_tbl = [emp_atr, sim_atr];
                site_tbl.Site = repmat(categorical(site_name), height(site_tbl), 1);
            
                self.tbl = [self.tbl; site_tbl]; 
            end
        
            formula = 'EmpAtrophy ~ SimAtrophy + (1|Site)';
            self.lme = fitlme(self.tbl, formula);
            self.set_lme();
        end
        
        function self = set_lme(self)
            predicted = fitted(self.lme, 'Conditional', false);
            empirical = self.lme.Variables.EmpAtrophy;
            self.lme_corr = corr(predicted, empirical);
            % disp(['correlation between observed and predicted atrophy: ', num2str(self.lme_corr)])

            lme_stats = anova(self.lme);
            self.pval = lme_stats.pValue(strcmp(lme_stats.Term, 'SimAtrophy'));
            % disp(['SimAtrophy pval: ', num2str(self.pval)])
        end

        %% Some plotting functions
        function plot_lme(self, params)
            figure; plotResiduals(self.lme, 'fitted');
            plot_fx(self, params);
            plot_predictions(self, params);
        end

        function plot_fx(self, params)
            n_sites = numel(params.site_names);
            predicted = fitted(self.lme);

            figure('Position', [400 250 1000 800]);
            for i = 1:n_sites
                site = string(self.tbl.Site);
                site_mask = site == string(params.site_names(i));
                subplot(2, 2, i);
                hold on;
                scatter(self.tbl.EmpAtrophy(site_mask), predicted(site_mask), 10);
                plot_bestfit(self.tbl.EmpAtrophy(site_mask), predicted(site_mask));
                hold off;
                xlabel('Empirical Atrophy');
                ylabel('Predicted Atrophy');
                title(params.site_names(i), ' Fixed Effects');
            end
        end

        function plot_predictions(self, params)
            predicted = fitted(self.lme, 'Conditional', false);
            predicted = reshape(predicted, params.n_rois, []);
            empirical = table2array(params.emp_atr);

            % Heatmap
            figure;
            subplot(1,2,1);
            imagesc(predicted);
            title('Predicted Atrophy');
            colorbar; clim([0 3.5]);
            subplot(1,2,2);
            imagesc(empirical);
            title('Empirical Atrophy');
            colorbar; clim([0 3.5]);
        end
    end
end