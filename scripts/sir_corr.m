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
    step_ = 1;

    % Compute correlation time series for all sites
    all_corrs_p = compute_corrs(params, sim_atrophy, step_, 'Pearson');

    % Find maximum correlations for each site
    gene = compute_max_corr(params, gene, all_corrs_p, step_);

    if params.null == "none"
        % Find time step where correlation reaches max
        gene.t_max = (find(all_corrs_p(:, 1) >= max(all_corrs_p(:, 1)), 1) - 1) * step_ + 1;
        % Find time step where correlation is greater than 0.5
        gene.t_05 = (find(all_corrs_p(:, 1) >= 0.5, 1) - 1) * step_ + 1;
    end

    if params.vis && params.null == "none"
        % Compute max atrophy per roi
        max_atr = compute_atr(params, sim_atrophy, gene.max_corr, step_);
        all_corrs_s = compute_corrs(params, sim_atrophy, step_, 'Spearman');

        % Show animated video plot
        % plot_video(params, sim_atrophy, all_corrs_p, all_corrs_s, step_, 2000);

        % Show static plots
        figure('Position', [100 100 500 1000], 'Color', 'white');
        subplot(2, 1, 1);
        plot_corrs(gene, all_corrs_p, all_corrs_s, step_);
        subplot(2, 1, 2);
        plot_scatter(params, gene, max_atr);
    end
end

function all_corrs = compute_corrs(params, sim_atrophy, step_, type)
    % returns (time series x sites) matrix
    % full time series only used for plot_corrs

    all_corrs = zeros(ceil(params.t/step_), params.n_sites);

    for i = 1:params.n_sites
        emp_atr = table2array(params.emp_atr(:, i));
        all_corrs(:, i) = corr(sim_atrophy(:, 1:step_:end), emp_atr, 'Type', type);
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

function plot_corrs(gene, all_corrs_p, all_corrs_s, step_)
    % produces correlation over the full time series
    corrs_p = all_corrs_p(:, 1);                       % lme beta only
    corrs_s = all_corrs_s(:, 1);                       % lme beta only
    % corrs = all_corrs;
 
    % t_total = length(corrs_p)*step_;
    t_total = 2000;  % limit x axis for vis
    hold on;
    plot(1:step_:t_total, corrs_p(1:t_total/step_), 'LineWidth', 1.5);
    plot(1:step_:t_total, corrs_s(1:t_total/step_), 'LineWidth', 1.5);
    % xline(gene.t_max, '-.r', num2str(corrs_p((gene.t_max-1)/step_)), 'LabelHorizontalAlignment', 'center', 'LabelVerticalAlignment','bottom');

    % t = title ({"Simulated and empirical atrophy", "risk gene: " +          ...
        % gene.risk_name, "clearance gene: " + gene.clear_name});
    t = title ({"Correlation over time"});
    % t.FontWeight = 'normal';
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
        scatter(sim_atr, emp_atr, 25, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.7);
        plot_bestfit(sim_atr, emp_atr, c=colors(i, :));
    end

    diagnosis_corr = num2str(table2array(gene.max_corr(1, 2)));
    % t = title({['Correlation = ', diagnosis_corr], ...
    %     "risk gene: " + gene.risk_name, "clearance gene: " + gene.clear_name});
    t = title ({"GMV at peak correlation"});
    % t.FontWeight = 'normal';
    xlabel("simulated GMV")
    ylabel('empirical GMV (z score)')
    hold off;
end

function plot_video(params, sim_atrophy, all_corrs_p, all_corrs_s, step_, max_t)
    % creates an animated video showing both scatter plot and correlation plot evolving through time
    % Left panel: correlation plot building up over time (Pearson and Spearman)
    % Right panel: scatter plot with moving dots as atrophy evolves
    % Automatically exports video to 'sir_simulation.mp4'
    
    % Setup figure
    fig = figure('Position', [100 100 500 500], 'Color', 'white');
    
    % Setup video writer
    v = VideoWriter('sir_simulation.avi', 'Motion JPEG AVI');
    v.FrameRate = 10;
    v.Quality = 90;
    open(v);
    
    % Calculate time points and data ranges
    n_timepoints = size(all_corrs_p, 1);
    
    max_timepoint = min(ceil(max_t/step_), n_timepoints);
    
    corrs_p = all_corrs_p(:, 1);  % Pearson, lme beta only
    corrs_s = all_corrs_s(:, 1);  % Spearman, lme beta only
    % emp_atr = params.emp_atr{:, 1};  % empirical atrophy
    
    % Pre-calculate scatter plot data only up to max_timepoint
    % sim_atr_data = zeros(params.n_rois, max_timepoint);
    % for t = 1:max_timepoint
    %     sim_atr_data(:, t) = normalize(sim_atrophy(:, (t-1)*step_+1));
    % end
    
    % Set up axis limits based on data range up to max_timepoint
    max_corr = max([max(corrs_p(1:max_timepoint)), max(corrs_s(1:max_timepoint))]);
    corr_ylim = [0, max(0.8, max_corr * 1.1)];
    % scatter_xlim = [min(sim_atr_data,[],1)' * 1.1, max(sim_atr_data,[],1)' * 1.1];
    % scatter_ylim = [min(emp_atr) * 1.1, max(emp_atr) * 1.1];

    colors = lines(3);            % match colors
    for t = 10:10:max_timepoint
        clf(fig);  % Clear figure
        
        % subplot: Correlation over time
        % subplot(2, 1, 1);
        hold on;
        
        % Plot correlations up to current time point
        time_axis = 1:step_:t*step_;
        plot(time_axis, corrs_p(1:t), 'Color', colors(1,:), 'LineWidth', 1.5, 'DisplayName', 'Pearson');
        plot(time_axis, corrs_s(1:t), 'Color', colors(2,:), 'LineWidth', 1.5, 'DisplayName', 'Spearman');

        % Highlight current points
        plot(t*step_, corrs_p(t), 'o', 'MarkerSize', 6, 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));
        plot(t*step_, corrs_s(t), 'o', 'MarkerSize', 6, 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', colors(2,:));

        legend({'Pearson', 'Spearman'}, 'Location', 'southeast');

        xlim([1, max_timepoint*step_]);
        ylim(corr_ylim);
        xlabel('t');
        ylabel('correlation');
        title(sprintf('Correlation over time (t = %d)', t*step_));
        hold off;
        
        % subplot: Scatter plot at current time
        % subplot(2, 1, 2);
        % hold on;
        % 
        % % Plot scatter
        % scatter(sim_atr_data(:, t), emp_atr, 25, colors(1,:), 'filled', 'MarkerFaceAlpha', 0.7);
        % plot_bestfit(sim_atr_data(:, t), emp_atr, c=colors(1, :));
        % current_corr = corr(sim_atr_data(:, t), emp_atr, 'Type', 'Pearson');
        % 
        % xlim(scatter_xlim(t, :));
        % ylim(scatter_ylim);
        % xlabel('simulated GMV');
        % ylabel('empirical GMV (z score)');
        % title(sprintf('GMV at time %d (r = %.2f)', t*step_, current_corr));
        % hold off;
        
        % pause(0.1); % Pause for visualization
        % Capture frame for video export
        frame = getframe(fig);
        writeVideo(v, frame);
        drawnow;
    end
    
    % Close video writer
    close(v);
end

