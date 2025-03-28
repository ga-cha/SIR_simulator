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
        % visualisation flag
        vis = false;

        %% input parameters (inside parenthesis are values used in Zheng 2019 paper)
        % model tuning parameters
        v = 1;                  % v: speed (1)
        dt = 0.01;              % dt: time step (0.01)
        t_total = 20000;        % T_total: total time steps (20000)
        init_number = 1;        % init_number: number of injected misfolded alpha-syn (1)
        prob_stay = 0.5;        % prob_stay: the probability of staying in the same region per unit time (0.5)
        trans_rate = 1;         % trans_rate: a scalar value, controlling the baseline infectivity
        beta_coeff;

        % parcellation specific parameters
        sconnLen;               % sconnLen: structural connectivity matrix (length) (estimated from HCP data)
        sconnDen;               % sconnDen: structural connectivity matrix (strength) (estimated from HCP data)
        ROIsize;                % ROIsize: region sizes (voxel counts)

        emp_atr = table();      % table indexed by site_name

        % derived parameters
        n_rois;                 % N_regions: number of regions (42)
        seed;                   % seed: seed region of misfolded alpha-syn injection (choose as you like? (^?^)= here substantia nigra)

        % optional null parameters
        null = "none";
        null_len;
        null_den;
        null_atr;
    end

    methods 
        function self = SIR_parameters(opt)
            if isfield(opt, 'vis'); self.vis = opt.vis; end
            self.null = opt.null;
            self.beta_coeff = opt.beta;

            self = self.set_params(opt);
            self = self.set_netw(opt);
            self = self.set_atrophy(opt);
        end

        function self = set_params(self, opt)
            ws = 'data/workspace_' + opt.parc + '.mat';
            load(ws, 'ROIsize');
            self.ROIsize = ROIsize;
            self.n_rois = length(ROIsize);
            if self.n_rois == 41         % DK + aseg
                self.seed = 40;
            elseif self.n_rois == 66     % Schaefer 100 + Tian S2
                self.seed = 51;
            elseif self.n_rois == 166    % Schaefer 300 + Tian S2
                self.seed = 151;
            else
                error("Accepts DK + aseg or Schaefer100/300 + Tian S2 only");
            end
        end

        function self = set_netw(self, opt)
            if opt.null ~= "rewired"
                ws = 'data/workspace_' + opt.parc + '.mat';
                load(ws, 'ifod_len_35', 'ifod_den_35');
                self.sconnLen = ifod_len_35;
                self.sconnDen = ifod_den_35;
            else
                ws = 'data/workspace_' + opt.parc + '_null.mat';
                load(ws, 'null_len', 'null_den');
                self.sconnLen = null_len;
                self.null_den = null_den;
            end
        end

        function self = set_atrophy(self, opt)
            if opt.null ~= "spatial"
                ws = 'data/workspace_' + opt.parc + '.mat';
                load(ws, 'emp_atr');
                self.emp_atr = emp_atr;
            else
                ws = 'data/workspace_' + opt.parc + '_null.mat';
                load(ws, 'null_atr');
                self.null_atr = null_atr;
                self.emp_atr = self.null_atr{1};
            end
        end

        function self = set_spatial(self, n)
            self.emp_atr = self.null_atr{n};
        end
    end
end