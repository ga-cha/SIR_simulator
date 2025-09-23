function out = calcResiduals(x, y, degree)
%CALCRESIDUALS calculate residuals of y wrt x

if nargin < 3; degree = 1; end

% p = polyfit(x, y, degree);
% yHat = polyval(p,x);
% out = y - yHat;

[p, ~, xfm] = polyfit(x, y, degree);
out = y - polyval(p, x, [], xfm);

end

