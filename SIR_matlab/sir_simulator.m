% SIRsimulator.m
%
% Package based off https://github.com/yingqiuz/SIR_simulator
% Zheng, Ying-Qiu, et al. PLoS biol. 17.11 (2019): e3000495.
%

function [Rnor_all, Rmis_all, Pnor_all, Pmis_all] = ...
    sir_simulator(params, gene, vis)

    % A function to simulate the spread of misfolded alpha-syn
    
    %% output parameters
    % Rnor_all: A N_regions * T_total matrix, recording the number of normal
    % alpha-syn in regions
    % Rmis_all: A N_regions * T_total matrix, recording the number of
    % misfolded alph-syn in regions
    % Pnor_all: a N_regions * N_regions * T_total matrix, recording the number of normal alpha-syn in paths
    % could be memory-consuming)
    % Pmis_all: a N_regions * N_regions * T_total matrix, recording the number of misfolded alpha-syn in paths
    % could be memory-consuming)
    
    % unpack parameters
    v = params.v;
    dt = params.dt;
    T_total = params.T_total;
    init_number = params.init_number;
    prob_stay = params.prob_stay;
    trans_rate = params.trans_rate;
    sconnDen = params.sconnDen;
    sconnLen = params.sconnLen;
    ROIsize = params.ROIsize;
    N_regions = params.N_regions;
    seed = params.seed;

    % set maximum timesteps for normal protein propagation
    % in Zheng 2019 this is different from T_total
    iter_max = 100000;
    
    
    % set the mobility pattern
    weights = sconnDen;
    
    % GC: A protein movement correction
    % the probability of exit from any given region is proportional to
    % edge contribution to total connectome weight
    prob_stay = 1 - rescale(sum(sconnDen)') .* prob_stay;
    
    weights = (1 - prob_stay) .* weights + prob_stay .* diag(sum(weights, 2)) ;
    
    % multinomial distribution
    % element (i,j) is the probability of moving from region i to edge (i,j)
    weights = weights ./ repmat(sum(weights, 2), 1, N_regions);
    weights(eye(N_regions, 'logical')) = 0;
    
    
    % convert gene expression (z) scores to probabilities
    clearance_rate = normcdf(zscore(gene.clear_gene));
    synthesis_rate = normcdf(zscore(gene.risk_gene));
    
    % store the number of normal/misfoled alpha-syn at each time step
    [Rnor_all, Rmis_all] = deal( zeros([N_regions, T_total]) );
    [Pnor_all, Pmis_all] = deal( zeros([N_regions, N_regions, T_total]) );
    % optionally store normal protein growth when filling network 
    if vis; [Rnor_nor_all] = deal(zeros([N_regions, iter_max])); end

    % Rnor, Rmis, Pnor, Pmis store results of single simulation at each time
    [Rnor, Rmis] = deal(zeros(N_regions, 1));   % number of normal/misfolded alpha-syn in regions
    [Pnor, Pmis] = deal(zeros(N_regions));      % number of normal/misfolded alpha-syn in paths
    
    % simplification of variables
    alphaTerm = (synthesis_rate .* ROIsize) .* dt;
    % betaTerm = exp(-(clearance_rate).*dt);
    % GC: A modification of clearance rate, since gene expression is a
    % proxy for protein expression
    betaTerm = exp(-((1/sqrt(2)) .* clearance_rate).*dt);

    sTerm = 1 ./ sconnLen .* dt .* v; sTerm(isinf(sTerm)) = 0;
    wTerm = weights .* dt;
    gamma0 = 1 .* trans_rate ./ ROIsize .* dt ; % the probability of getting misfolded
    
    
    %% normal alpha-syn growth
    % fill the network with normal proteins
    % typically converges in t < 100000
    % disp('normal alpha synuclein growth');
    for t = 1:iter_max
        %%% moving process
        % regions towards paths
        % movDrt stores the number of proteins towards each region. i.e.
        % element in kth row lth col denotes the number of proteins in region k
        % moving towards l
        movDrt = Rnor .* wTerm; % implicit expansion
    
        % paths towards regions
        % update moving
        movOut = Pnor .* sTerm; % longer path & smaller v = lower probability of moving out of paths
    
        Pnor = Pnor - movOut + movDrt;
        Rtmp = Rnor;
        Rnor = Rnor + sum(movOut, 1)' - sum(movDrt, 2);
    
        %%% growth process
        Rnor = Rnor .* betaTerm + alphaTerm;
    
        if vis; Rnor_nor_all(:, t) = Rnor; end
    
        if abs(Rnor - Rtmp) < (1e-7 * Rtmp); break; end
    end
    %% misfolded protein spreading process
    % inject misfolded alpha-syn
    Rmis(seed) = init_number;
    % disp('misfolded alpha synuclein spreading');
    for t = 1:T_total
        %%% moving process
        % normal proteins: region -->> paths
        movDrt_nor = Rnor .* wTerm; % implicit expansion
    
        % normal proteins: paths -->> regions
        movOut_nor = Pnor .* sTerm;
    
        % misfolded proteins: region -->> paths
        movDrt_mis = Rmis .* wTerm; % implicit expansion
    
        % misfolded proteins: paths -->> regions
        movOut_mis = Pmis .* sTerm;
    
        % update regions and paths
        Pnor = Pnor - movOut_nor + movDrt_nor;
        Rnor = Rnor + sum(movOut_nor, 1)' - sum(movDrt_nor, 2);
    
        Pmis = Pmis - movOut_mis + movDrt_mis;
        Rmis = Rmis + sum(movOut_mis, 1)' - sum(movDrt_mis, 2);
    
        misProb = 1 - exp( -Rmis .* gamma0 ) ; % trans_rate: default
        % number of newly infected
        N_misfolded = Rnor .* exp(-clearance_rate) .* misProb ;
        
        % update
        Rnor = Rnor .* betaTerm + alphaTerm - N_misfolded;
        Rmis = Rmis .* betaTerm             + N_misfolded;
    
        Rnor_all(:, t) = Rnor ;
        Rmis_all(:, t) = Rmis ;
    
        % uncomment the following lines if you want outputs of alpha-syn in
        % paths
        % Pnor_ave(:, :, t) = Pnor;
        % Pmis_ave(:, :, t) = Pmis;
    end

    if vis; plot_propagation(Rnor_nor_all, Rnor_all, Rmis_all); end

end

function plot_propagation (Rnor_nor_all, Rnor_all, Rmis_all)
    figure;
    plot(Rnor_nor_all');
    t = title ("Normal protein propagation");
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("Ni");

    figure;
    plot(Rnor_all');
    t = title ("Misfolded protein propagation: normal");
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("Ni");

    figure;
    plot(Rmis_all');
    t = title ("Misfolded protein propagation: misfolded");
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("Mi");
    ylim([0 30000])
end