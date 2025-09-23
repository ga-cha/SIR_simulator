% clear_names = "RNF144A";
% opt.risk_names = ["TCAP"];
% opt.parc = "S132" ;
% opt.init = 1;
% opt.dbg = false;
% opt.vis = false;
% opt.pf = false;
% opt.null = "none";

%% 
% tic

n = 1e-13;
t = zeros(66, 1);
for i = 1:66
    for j = 1:20
        init = n*sqrt(10)^j;
        gene_corrs = main("RNF144A", risk_names=["TCAP"], parc ="S132", dbg=true, seed=i, init=init, pf=false);
        if (~isnan(gene_corrs{:, "correlation"}))
            break
        end
    end
    t(i) = init;
end

% toc
% 
% figure('Color', 'white')
% plot(t, '.', 'LineStyle', 'none')
% set(gca, 'YScale', 'log')
% t = title ("Spread threshold (RNF144A/TCAP)");
% t.FontWeight = 'normal';
% xlabel("seed region");
% ylabel("misfolded seed quantity");
