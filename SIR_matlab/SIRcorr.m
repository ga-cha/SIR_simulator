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

function [bgs_max, cobre_max, hcpep_max, stages_max, tstep, err] = ...
    SIRcorr(params, gene, n)

    % Take log2fold change 
    % gene.sim_atrophy = real(log2(gene.sim_atrophy));
    if params.null == "spatial"; params = set_spatial(params, n); end
    % Calculate position and value of max correlation coefficient
    bgs_corrs = corr(gene.sim_atrophy, params.bgs, 'Type', 'Pearson');
    cobre_corrs = corr(gene.sim_atrophy, params.cobre, 'Type', 'Pearson');
    hcpep_corrs = corr(gene.sim_atrophy, params.hcpep, 'Type', 'Pearson');
    stages_corrs = corr(gene.sim_atrophy, params.stages, 'Type', 'Pearson');
    % bgs_corrs = corr(gene.sim_atrophy, params.bgs, 'Type', 'Spearman');
    % cobre_corrs = corr(gene.sim_atrophy, params.cobre, 'Type', 'Spearman');
    % hcpep_corrs = corr(gene.sim_atrophy, params.hcpep, 'Type', 'Spearman');
    % stages_corrs = corr(gene.sim_atrophy, params.stages, 'Type', 'Spearman');
    % Truncate first 1000 timepoints which heavily weight seed region
    [bgs_max, tstep] = max(bgs_corrs(1001:end));
    cobre_max = max(cobre_corrs(1001:end));
    hcpep_max = max(hcpep_corrs(1001:end));
    stages_max = max(stages_corrs(1001:end));
    tstep = tstep + 1000;
    err = 0;
    % err= calcResiduals(sim_atrophy(:, tstep), params.bgs);

    if params.vis
        plot_corrs(bgs_corrs, cobre_corrs, hcpep_corrs, gene);
        max_corr = (bgs_max + cobre_max + hcpep_max)./3;
        plot_scatter(tstep, params, gene, max_corr); 
    end
end

function plot_corrs(bgs_corrs, cobre_corrs, hcpep_corrs, gene)
    max_t = 10000;
    figure;
    hold on;
    plot([bgs_corrs(1:max_t), cobre_corrs(1:max_t), hcpep_corrs(1:max_t)]);
    t = title ({"Simulated and empirical atrophy", "risk gene: " + gene.risk_name, "clearance gene: " + gene.clear_name});
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("correlation");
    legend("BGS", "COBRE", "HCPEP");
end

function plot_scatter(tstep, params, gene, max_corr)
    emp_atrophy = (params.bgs + params.cobre + params.hcpep)./3;
    figure;
    hold on;
    scatter(gene.sim_atrophy(:, tstep), emp_atrophy); 

    % polyfit fits a degree 1 polynomial
    % polyval evaluates polynomial at specified points
    [p, ~, xfm] = polyfit(gene.sim_atrophy(:, tstep), emp_atrophy, 1);
    y = polyval(p, gene.sim_atrophy(:, tstep), [], xfm);
    plot(gene.sim_atrophy(:, tstep), y');

    t = title({['Pearson correlation = ', num2str(max_corr)], "risk gene: " + gene.risk_name, ...
        "clearance gene: " + gene.clear_name});
    % t = title({['Spearman correlation = ', num2str(max_corr)], "risk gene: " + gene.risk_name, ...
    %     "clearance gene: " + gene.clear_name});
    t.FontWeight = 'normal';
    xlabel("simulated atrophy")
    ylabel('empirical atrophy (z score)')
end