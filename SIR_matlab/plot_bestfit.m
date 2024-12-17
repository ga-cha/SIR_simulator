%
% plot_bestfit
%
% plots a line of best fit
%
%

function plot_bestfit(x, y, opt)
    arguments
        x;
        y;
        opt.c = lines(1);
        opt.d = 1;
    end

    % assert (nargin > 2, "Usage: plot_best_fit(x, y, c=color, " + ...
    %     "d=polynomial degree");


    % polyfit fits a degree 1 polynomial
    % polyval evaluates polynomial at specified points
    [p, ~, xfm] = polyfit(x, y, opt.d);
    y2 = polyval(p, x, [], xfm);
    plot(x, y2', 'color', opt.c);
end
