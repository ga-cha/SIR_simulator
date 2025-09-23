% SIRparameters.m
%
% Gabriella Chan 25/06/24
% gabriella.chan@monash.edu
%
% Introduce Parameter Object for SIR simulator main
% Holds model tuning parameters, parcellation data and atrophy data
% 

classdef SIR_parameters
    properties
        %% input parameters (inside parenthesis are values used in Zheng 2019 paper)
        % model tuning parameters
        v = 1;                  % v: speed (1)
        dt = 0.01;              % dt: time step (0.01)
        t = 10000;              % T_total: total time steps (20000)
        init_number;            % init_number: number of injected misfolded alpha-syn (1)
        prob_exit = 0.1;        % prob_stay: the probability of staying in the same region per unit time (0.5)
        trans_rate = 1;         % trans_rate: a scalar value, controlling the baseline infectivity
        beta_coeff = 1/sqrt(2);

        % user input parameters
        sc_length;              % sconnLen: structural connectivity matrix (length) (estimated from HCP data)
        sc_weight;              % sconnDen: structural connectivity matrix (strength) (estimated from HCP data)
        roi_size;               % ROIsize: region sizes (voxel counts)
        n_rois;                 % N_regions: number of regions (42)
        seed;                   % seed: seed region of misfolded alpha-syn injection (choose as you like? (^?^)= here substantia nigra)
        emp_atr = table();      % empirical atrophy, indexed by site_name
        n_sites;                % derived from emp_atr, here for convenience
        site_names;             % derived from emp_atr, here for convenience

        % flags and null parameters
        vis;                    % visualisation flag
        pf;                     % parfor flag
        null = "none";          % null flag
        null_weight;            % rewired null sc weights
        null_atr;               % spatial null empirical atrophy

        ws;
    end

    methods 
        function self = SIR_parameters(opt)
            self.ws = '../data/workspace_' + opt.parc + '.mat';

            self.vis = opt.vis;
            self.pf = (opt.pf && ~opt.vis); % visualisation unavailable with parfor
            self.null = opt.null;
                    
            load(self.ws, 'roi_size');
            self.roi_size = roi_size;
            self.n_rois = length(roi_size);
            self.init_number = opt.init;

            self = self.set_seed(opt);
            self = self.set_netw(self.ws);
            self = self.set_emp_atr(self.ws);
            
            self.n_sites = width(self.emp_atr);
            self.site_names = string(self.emp_atr.Properties.VariableNames);
        end

        function self = set_seed(self, opt)
            % Sets seed, default left hippocampus
            if isfield(opt, 'seed')
                assert(opt.seed <= self.n_rois, ...
                    "Seed is not set to a valid region");
                self.seed = opt.seed;
            else
                if self.n_rois == 41         % DK + aseg
                    self.seed = 41;
                elseif self.n_rois == 66     % Schaefer 100 + Tian S2
                    self.seed = 51;
                elseif self.n_rois == 166    % Schaefer 300 + Tian S2
                    self.seed = 151;
                else
                    error("Accepts DK + aseg or Schaefer100/300 + Tian S2 only");
                end
            end
        end

        function self = set_netw(self, ws)
            load(ws, 'ifod_l_35', 'ifod_w_35');
            self.sc_length = ifod_l_35;
            self.sc_weight = ifod_w_35;
        end

        function self = set_emp_atr(self, ws)
            load(ws, 'emp_atr');
            self.emp_atr = emp_atr;
        end

        function self = set_null(self)
            if self.null == "spatial"
                load(self.ws, 'null_atr');
                self.null_atr = null_atr;
                self.emp_atr = self.null_atr{1};
            elseif self.null == "rewired"
                load(self.ws, 'null_l', 'null_w');
                self.sc_length = null_l;
                self.null_weight = null_w;
                self.emp_atr = self.emp_atr(:,1);
            end
            self.n_sites = width(self.emp_atr);
            self.site_names = string(self.emp_atr.Properties.VariableNames);
        end

        function self = set_spatial(self, n)
            self.emp_atr = self.null_atr{n};
        end
    end
end