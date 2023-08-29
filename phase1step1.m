function [Rnor, Rtmp, Pnor, Ptmp, movDrt, movOut] = phase1step1(Rnor, Pnor, alphaTerm, betaTerm, sTerm, wTerm, n)

for ii = 1:n 
    movDrt = Rnor .* wTerm; % IMPLICIT EXPANSION

    % paths towards regions
    % update moving
    movOut = Pnor .* sTerm; % longer path & smaller v = lower probability of moving out of paths

    Ptmp = Pnor;
    Pnor = Pnor - movOut + movDrt;
    Rtmp = Rnor;
    Rnor = Rnor + sum(movOut, 1)' - sum(movDrt, 2);

    %%% growth process
    Rnor = Rnor.*betaTerm + alphaTerm; % corresponds to Rtest2
%     Rnor = Rnor.*betaTerm + alphaTerm + sum(movOut, 1)' - sum(movDrt, 2); % corresponds to Rtest ;

end

end
