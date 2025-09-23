% SIRcorr
%
% Gabriella Chan 28/06/23
% gabriella.chan@monash.edu
% Monash University
%
% Given simulated atrophy, this module calculates similarity with
% empirical atrophy.
% Simulated atrophy input contains atrophy per region per timepoint. For
% each timepoint, then, we can calculate empirical correlation. We find
% the maximum and optionally plot

function gene = sir_corr(params, gene, sim_atrophy)

    % evaluate at step_ sized intervals
    % TODO: double check for off-by-ones
    step_ = 3;

    % Compute correlation time series for all sites
    all_corrs = compute_corrs(params, sim_atrophy, step_);

    % Find maximum correlations for each site
    gene = compute_max_corr(params, gene, all_corrs, step_);

    if params.null == "none"
        % Find time step where correlation reaches max
        gene.t_max = (find(all_corrs(:, 1) >= max(all_corrs(:, 1)), 1) - 1) * step_ + 1;
        % Find time step where correlation is greater than 0.5
        gene.t_05 = (find(all_corrs(:, 1) >= 0.5, 1) - 1) * step_ + 1;
    end

    if params.vis && params.null == "none"
        % Compute max atrophy per roi
        max_atr = compute_atr(params, sim_atrophy, gene.max_corr, step_);
        % Plot visualisation
        % figure('Position', [100 100 1500 600], 'Color', 'white');
        figure('Position', [100 100 1000 400], 'Color', 'white');
        subplot(1, 2, 1);
        plot_corrs(gene, all_corrs, step_);
        subplot(1, 2, 2);
        plot_scatter(params, gene, max_atr);
    end
end

function all_corrs = compute_corrs(params, sim_atrophy, step_)
    % returns (time series x sites) matrix
    % full time series only used for plot_corrs

    all_corrs = zeros(ceil(params.t/step_), params.n_sites);

    for i = 1:params.n_sites
        emp_atr = table2array(params.emp_atr(:, i));
        all_corrs(:, i) = corr(sim_atrophy(:, 1:step_:end), emp_atr, 'Type', 'Pearson');
    end
end

function gene = compute_max_corr(params, gene, all_corrs, step_)
    % populates gene.max_corr ([site_name, site_max, site_t] x sites) table
    % takes full time series across sites and find maximum correlation

    % skip_ = ceil(1000/step_);       % Ignore first 1000 heavily weighted seed regions
    skip_ = 0;

    for i = 1:params.n_sites
        [site_max, site_t] = max(all_corrs(skip_+1:end, i)); 
        gene.max_corr(i, :) = {params.site_names(i), site_max, site_t + skip_};
    end
end

function max_atr = compute_atr(params, sim_atrophy, max_corr, step_)
    % returns (rois x sites) matrix of simulated atrophy at peak fit
    % only used in plot_scatter
    max_atr = zeros(params.n_rois, params.n_sites);

    for i = 1:params.n_sites
        max_atr(:, i) = sim_atrophy(:, (max_corr.site_t(i)-1)*step_+1);
    end
end

function plot_corrs(gene, all_corrs, step_)
    % produces correlation over the full time series
    corrs = all_corrs(:, 1);                       % lme beta only
    % corrs = all_corrs;
    
    plot(1:3:length(corrs)*step_, corrs);

    xline(gene.t_max, '-.r', num2str(corrs((gene.t_max-1)/step_)), 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','bottom');

    % t = title ({"Simulated and empirical atrophy", "risk gene: " +          ...
        % gene.risk_name, "clearance gene: " + gene.clear_name});
    t = title ({"Simulation correlation over time"});
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("correlation");
    ylim([0 0.8]);
end

function plot_scatter(params, gene, max_atr)
    % produces a scatter plot and line of best fit
    % between simulated and empirical atrophy
    hold on;
    n_sites = 1;                                     % lme beta only
    % n_sites = params.n_sites;                      % lme beta and all sites
    colors = lines(n_sites);                         % match color of scatter and best fit
    for i = 1:n_sites
        sim_atr = normalize(max_atr(:, i));          % rescale for vis
        emp_atr = params.emp_atr{:, i};
        scatter(sim_atr, emp_atr, 25, colors(i,:));
        plot_bestfit(sim_atr, emp_atr, c=colors(i, :));
    end

    diagnosis_corr = num2str(table2array(gene.max_corr(1, 2)));
    % t = title({['Correlation = ', diagnosis_corr], ...
    %     "risk gene: " + gene.risk_name, "clearance gene: " + gene.clear_name});
    t = title ({"GMV at peak correlation"});
    t.FontWeight = 'normal';
    xlabel("simulated GMV")
    ylabel('empirical GMV (z score)')
    hold off;
end