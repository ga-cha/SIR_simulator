% main.m
%                                                                                                                                    
%
% Package based off https://github.com/yingqiuz/SIR_simulator
% Zheng, Ying-Qiu, et al. PLoS biol. 17.11 (2019): e3000495.
%
% Gabriella Chan 30/06/23
% gabriella.chan@monash.edu
% Monash University
%
% a script to simulate atrophy accrual due to the accumulation of misfolded
% alpha-syn aggregates
% 42-region parcellation
% load gene expressions, real atrophy, ROIsize, functional connectivity...
%
%

% load('data_gc/42regions/GC_workspace.mat');
%load('data_yqz/42regions/YQZ_workspace.mat');
% load structural connectivity
%load('data_yqz/42regions/sc35.mat');

N_regions = 41;
v = 1;
dt = 0.01;
T_total = 20000;
init_number = 1;
syn_control = ROIsize;
prob_stay = 0.5;
trans_rate = 1;
% init seed to hip
seed = 40;
% load your GBA, SNCA, sconnDen, sconnLen, ROISize ....

[gene_corrs, sim_atrophy] = SIRiterator(N_regions, v, dt, T_total, parvalb, genes, sconnLen, ...
    sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, ...
    trans_rate, emp_atrophy);
tail(gene_corrs)
writetable(gene_corrs, "gene_corrs_parvalb_clear.csv")