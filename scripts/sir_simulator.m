% SIRsimulator.m
%
% Package based off https://github.com/yingqiuz/SIR_simulator
% Zheng, Ying-Qiu, et al. PLoS biol. 17.11 (2019): e3000495.
%

function proteins = sir_simulator(params, gene, vis)

    % A function to simulate the spread of misfolded alpha-syn
    % unpack parameters
    v = params.v;
    dt = params.dt;
    t = params.t;
    init_number = params.init_number;
    prob_exit = params.prob_exit;
    trans_rate = params.trans_rate;
    sc_weight = params.sc_weight;
    sc_length = params.sc_length;
    roi_size = params.roi_size;
    n_rois = params.n_rois;
    seed = params.seed;
    beta_coeff = params.beta_coeff;
    
    % set maximum timesteps for normal protein propagation
    iter_max = 1000000;         % (1000000000), typically converges in t < 100000
    
    %% Precompute mobility and synthesis/clearance rates
    % GC: Adjust probabilities based on structural connectivity
    scaled_sc = sum(sc_weight)' ./ max(sum(sc_weight));
    prob_stay = 1 - scaled_sc .* prob_exit;

    % Compute transition probabilities
    weights = (1 - prob_stay) .* sc_weight + prob_stay .* diag(sum(sc_weight, 2));
    weights = weights ./ sum(weights, 2);       % Normalize rows
    weights(eye(n_rois, 'logical')) = 0;        % Remove self-loops
 
    % Convert gene expression scores to probabilities
    clearance_rate = normcdf(zscore(gene.clear_gene));
    synthesis_rate = normcdf(zscore(gene.risk_gene));

    % Precompute constants for growth and movement
    alpha_ = (synthesis_rate .* roi_size) .* dt;        % Synthesis term
    beta_ = exp(-(beta_coeff .* clearance_rate) .* dt); % Clearance term
    s_ = 1 ./ sc_length .* dt .* v; s_(isinf(s_)) = 0;  % Path movement term
    w_ = weights .* dt;                                 % Region-to-path movement term
    gamma_ = trans_rate ./ roi_size .* dt;              % Misfolding probability

    % Initialize storage for results
    r_fill = zeros(n_rois, iter_max);           % Steady-state normal protein growth
    [r_nor, r_mis] = deal(zeros(n_rois, t));    % Normal/misfolded proteins over time
    [e_nor, e_mis] = deal(zeros(n_rois));       % Normal/misfolded proteins in paths

    %% Simulate normal protein growth
    r_fill(:, 1) = alpha_;                      % Initialize the first column of r_fill
    for i = 2:iter_max
        % Movement: regions -> paths -> regions
        movDrt = r_fill(:, i-1) .* w_;          % Region-to-edge 
        movOut = e_nor .* s_;                   % Edge-to-region 
        e_nor = e_nor - movOut + movDrt;
        r_fill(:, i) = r_fill(:, i-1) + sum(movOut, 1)' - sum(movDrt, 2);

        % Growth
        r_fill(:, i) = r_fill(:, i) .* beta_ + alpha_;

        % Check for convergence
        if all(abs(r_fill(:, i) - r_fill(:, i-1)) < (1e-7 * r_fill(:, i-1)))
            r_nor(:, 1) = r_fill(:, i);         % Set initial conditions
            break;
        end
    end

    %% Simulate misfolded protein spreading
    r_mis(seed, 1) = init_number;  % Inject misfolded proteins into the seed region
    for i = 1:t-1
        % Movement: normal proteins
        movDrt_nor = r_nor(:, i) .* w_;         % Region-to-edge 
        movOut_nor = e_nor .* s_;               % Edge-to-region 

        % Movement: misfolded proteins
        movDrt_mis = r_mis(:, i) .* w_;         % Region-to-edge 
        movOut_mis = e_mis .* s_;               % Edge-to-region 

        % Update regions and paths
        e_nor = e_nor - movOut_nor + movDrt_nor;
        r_nor(:, i+1) = r_nor(:, i) + sum(movOut_nor, 1)' - sum(movDrt_nor, 2);

        e_mis = e_mis - movOut_mis + movDrt_mis;
        r_mis(:, i+1) = r_mis(:, i) + sum(movOut_mis, 1)' - sum(movDrt_mis, 2);

        % Misfolding process
        mis_prob = 1 - exp(-r_mis(:, i) .* gamma_);  % Probability of misfolding
        n_mis = r_nor(:, i) .* exp(-clearance_rate) .* mis_prob;  % Newly misfolded proteins

        % Update states
        r_nor(:, i+1) = r_nor(:, i+1) .* beta_ + alpha_ - n_mis;
        r_mis(:, i+1) = r_mis(:, i+1) .* beta_ + n_mis;
    end

    proteins.r_nor = r_nor;
    proteins.r_mis = r_mis;

    if vis && params.null == "none"
        plot_propagation(r_fill, r_nor, r_mis);
    end

end

function plot_propagation (r_fill, r_nor, r_mis)

    figure('Position', [200 200 1500 600], 'Color', 'white');
    subplot(1, 3, 1);
    plot(r_fill');
    t = title ("Steady state");
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("N_i");

    subplot(1, 3, 2);
    plot(r_nor');
    t = title ("N_i");
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("N_i");

    subplot(1, 3, 3);
    plot(r_mis');
    t = title ("M_i");
    t.FontWeight = 'normal';
    xlabel("t");
    ylabel("M_i");
    ylim([0 30000])
end


    %% New approach to mapping steady state: 
    % regions and paths both become entities in a markov chain, with
    % distinct transition probabilities

    % tic; 
    % 
    % sconnMask = logical(sconnDen); 
    % [edgeX, edgeY] = find(sconnMask);
    % N_edges = height(edgeX); 
    % 
    % mat1 = sparse(1:N_regions, 1:N_regions, prob_stay); 
    % 
    % mat2 = sparse(edgeX, 1:N_edges, 1./sconnLen(sconnMask));
    % 
    % sconnDenStr = sum(sconnDen)'; 
    % mat3 = sparse((1:N_edges)', edgeX, ...
    %     (1-prob_stay(edgeX)) .* sconnDen(sconnMask) ./ sconnDenStr(edgeX) ); 
    % 
    % mat4 = sparse(1:N_edges, 1:N_edges, 1 - 1./sconnLen(sconnMask)); 
    % 
    % mat = [mat1, mat2; mat3, mat4];
    % 
    % % prev = zeros(N_regions + N_edges, 1); 
    % % for t = 1:iter_max
    % %     next = mat * prev; 
    % %     next(1:N_regions) = next(1:N_regions) .* betaTerm + alphaTerm; 
    % % 
    % %     if abs(next - prev) < (1e-7 * prev); break; end
    % % 
    % %     prev = next;
    % % end
    % [next,~] = eigs(mat, 1, 1); 
    % 
    % Rnor = next(1:N_regions); 
    % Pnor = full(sparse( edgeX, edgeY, next(N_regions+1:end) ));
    % 
    % toc; 
    % 
    % [Rnor, Rmis] = deal(zeros(N_regions, 1));   % number of normal/misfolded alpha-syn in regions
    % [Pnor, Pmis] = deal(zeros(N_regions));      % number of normal/misfolded alpha-syn in paths
    
    