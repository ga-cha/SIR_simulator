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

function [bgs_max, cobre_max, hcpep_max, stages_max, tstep, err] = SIRcorr(sim_atrophy, emp_atrophy_bgs, ...
    emp_atrophy_cobre, emp_atrophy_hcpep, emp_atrophy_stages, risk_name, clear_name, plotting)

    % Take log2fold change 
    sim_atrophy = real(log2(sim_atrophy));

    % Calculate position and value of max correlation coefficient
    [bgs_max, tstep] = corry(sim_atrophy, emp_atrophy_bgs, plotting, risk_name, clear_name);
    [cobre_max, ~] = corry(sim_atrophy, emp_atrophy_cobre, false, risk_name, clear_name);
    [hcpep_max, ~] = corry(sim_atrophy, emp_atrophy_hcpep, false, risk_name, clear_name);
    [stages_max, ~] = corry(sim_atrophy, emp_atrophy_stages, false, risk_name, clear_name);

    err= calcResiduals(sim_atrophy(:, tstep), emp_atrophy_bgs);
    % err = calcResiduals(sim_atrophy(:, tstep), emp_atrophy);

    if plotting
        figure;
        hold on;
        scatter(sim_atrophy(:, tstep), emp_atrophy_bgs); 
    
        % polyfit fits a degree 1 polynomial
        % polyval evaluates polynomial at specified points
        p = polyval(polyfit(sim_atrophy(:, tstep), emp_atrophy_bgs, 1), sim_atrophy(:, tstep));
        plot(sim_atrophy(:, tstep), p);

        t = title(['Spearman correlation = ', num2str(bgs_max)]);
        % t = title(['Pearson correlation = ', num2str(bgs_max)]);
        t.FontWeight = 'normal';
        xlabel({"simulated atrophy", "risk gene: " + risk_name, "clearance gene: " + clear_name})
        ylabel('empirical atrophy (z score)')
    end
end

function [max_corr, tstep] = corry(sim_atrophy, emp_atrophy, plotting, risk_name, clear_name)
    corrs = corr(sim_atrophy, emp_atrophy, 'Type', 'Spearman');
    % corrs = corr(sim_atrophy, emp_atrophy);
    if plotting
        figure;
        plot(corrs);
        t = title ({"risk gene: " + risk_name, "clearance gene: " + clear_name});
        t.FontWeight = 'normal';
        xlabel("t")
        ylabel("correlation")
    end
    % The first 1000 timepoints before spreading heavily weight the seed
    % region. Truncated for checking 
    [max_corr, tstep] = max(corrs(1001:end));
    tstep = tstep + 1000;
end
