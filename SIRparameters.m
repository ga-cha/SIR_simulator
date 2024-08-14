% SIRparameters.m
%
% Gabriella Chan 25/06/24
% gabriella.chan@monash.edu
%
% Introduce Parameter Object for SIR simulator main
% Holds model tuning parameters, parcellation data and atrophy data
% 

classdef SIRparameters
    properties
        vis = false;

        %% input parameters (inside parenthesis are values used in Zheng 2019 paper)
        % model tuning parameters
        v = 1;                  % v: speed (1)
        dt = 0.01;              % dt: time step (0.01)
        T_total = 20000;        % T_total: total time steps (20000)
        init_number = 1;        % init_number: number of injected misfolded alpha-syn (1)
        prob_stay = 0.5;        % prob_stay: the probability of staying in the same region per unit time (0.5)
        trans_rate = 1;         % trans_rate: a scalar value, controlling the baseline infectivity

        % parcellation specific parameters
        sconnLen;               % sconnLen: structural connectivity matrix (length) (estimated from HCP data)
        sconnDen;               % sconnDen: structural connectivity matrix (strength) (estimated from HCP data)
        ROIsize;                % ROIsize: region sizes (voxel counts)
        bgs;
        cobre;
        hcpep;
        stages;

        % derived parameters
        N_regions;              % N_regions: number of regions (42)
        seed;                   % seed: seed region of misfolded alpha-syn injection (choose as you like? (^?^)= here substantia nigra)

    end

    methods 
        function self = set_atrophy(self, bgs, cobre, hcpep, stages)
            self.bgs = bgs;
            self.cobre = cobre;
            self.hcpep = hcpep;
            self.stages = stages;
        end

        function self = set_netw(self, sconnLen, sconnDen, ROIsize)
            self.sconnLen = sconnLen;
            self.sconnDen = sconnDen;
            self.ROIsize = ROIsize;

            self.N_regions = length(ROIsize);

            if self.N_regions == 41         % DK + aseg
                self.seed = 40;
            elseif self.N_regions == 66     % Schaefer 100 + Tian S2
                self.seed = 51;
            elseif self.N_regions == 166    % Schaefer 300 + Tian S2
                self.seed = 151;
            else
                error("Accepts DK + aseg or Schaefer100/300 + Tian S2 only");
            end           
        end
    end
end