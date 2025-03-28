% SIRcorr
%
% Gabriella Chan 28/06/23
% gabriella.chan@monash.edu
% Monash University
%
% Given simulated atrophy, this module calculates similarity with
% empirical atrophy.
% Simulated atrophy input contains atrophy per region per timepoint. For
% each timepoint, then, we can calculate empirical correlation. We plot
%  and return the maximum

function gene = sir_corr(params, gene, sim_atrophy)

    % Take log2fold change 
    % sim_atrophy = real(log2(sim_atrophy));

    gene.max_atr = table();
    gene.max_corr = table();
    all_corrs = table();

    for i = 1:width(params.emp_atr)
        site_name = string(params.emp_atr.Properties.VariableNames(i));
        emp_atr = table2array(params.emp_atr(:, i));
        site_corrs = corr(sim_atrophy, emp_atr, 'Type','Pearson');
        [site_max, site_t] = max(site_corrs(1001:end));     % first 1000 heavily weight seed region
        
        % full time series of correlations (time series x sites matrix)
        % TODO: this is only used for plotting, move it down(?)
        site_corrs = table(site_corrs, 'VariableNames', site_name);
        all_corrs = [all_corrs, site_corrs];                                %#ok<AGROW>

        % peak correlation i.e. | Site Name | Max Correlation | Time Step |
        site_corr = table(site_name, site_max, site_t+1000, ...
            'VariableNames', {'SiteName', 'SiteCorr', 'SiteT'});
        gene.max_corr = [gene.max_corr; site_corr];

        % simulated atrophy at peak correlation (rois x sites matrix)
        site_atr = table(sim_atrophy(:, site_t+1000), 'VariableNames', site_name);
        gene.max_atr = [gene.max_atr, site_atr];
    end

    if params.vis
        figure('Position', [10 10 1500 600], 'Color', 'white');
        subplot(1, 2, 1);
        plot_corrs(gene, all_corrs);
        subplot(1, 2, 2);
        plot_scatter(gene, params); 
    end
end
    
function plot_corrs(gene, all_corrs)
    % corrs = table2array(all_corrs);
    % corrs = table2array(all_corrs(:, 1));
    corrs = [table2array(all_corrs(:, 1)), table2array(all_corrs(:, 22:end))];
    
    plot(corrs);
    t = title ({"Simulated and empirical atrophy", "risk gene: " +          ...
        gene.risk_name, "clearance gene: " + gene.clear_name});
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("correlation");
    % legend(all_corrs.Properties.VariableNames(:));
end

function plot_scatter(gene, params)
    hold on;
    % n_sites = width(gene.max_atr);
    n_sites = 1;
    colors = lines(n_sites);    % match color of scatter and best fit
    for i = 1:n_sites
        sim_atr = gene.max_atr{:, i};
        sim_atr = normalize(sim_atr);       % rescale for vis
        emp_atr = params.emp_atr{:, i};
        scatter(sim_atr, emp_atr, 25, colors(i,:), ...
            'DisplayName', gene.max_atr.Properties.VariableNames{i});
        plot_bestfit(sim_atr, emp_atr, c=colors(i, :));
    end

    diagnosis_corr = num2str(table2array(gene.max_corr(1, 2)));
    t = title({['Correlation = ', diagnosis_corr], ...
        "risk gene: " + gene.risk_name, "clearance gene: " + gene.clear_name});
    t.FontWeight = 'normal';
    xlabel("simulated atrophy")
    ylabel('empirical atrophy (z score)')
    hold off;
end

% figure;
% hold on;
% Site = self.tbl.Site;
% colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];
% gscatter(self.tbl.EmpAtrophy, self.tbl.SimAtrophy, self.tbl.Site);
% for i = 1:4
%     site_name = site_names{i};
%     emp_atrophy = self.tbl.EmpAtrophy(self.tbl.Site == site_name);
%     sim_atrophy = self.tbl.SimAtrophy(self.tbl.Site == site_name);
%     [p, ~, xfm] = polyfit(emp_atrophy, sim_atrophy, 1);
%     y = polyval(p, emp_atrophy, [], xfm);
%     plot(emp_atrophy, y', 'color', colors(i));
% end
% hold off;