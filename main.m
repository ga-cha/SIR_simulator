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

load('data/42regions/GC_workspace.mat');
load('data/42regions/some_genes.mat');
% load structural connectivity
load('data/42regions/sc35.mat');

N_regions = 42;
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

gene_corrs = SIRiterator(N_regions, v, dt, T_total, GBA, genes, sconnLen, ...
    sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, ...
    trans_rate, emp_atrophy);
head(gene_corrs)

