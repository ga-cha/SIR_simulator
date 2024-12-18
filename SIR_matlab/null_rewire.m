% null_rewire
%
% Gabriella Chan 01/08/23
% gabriella.chan@monash.edu
% Monash University
%
% This module calls randmio_und.m from the Brain Connectivity Toolbox to
% generate rewired null networks using the Maslov-Sneppen algorithm
% https://sites.google.com/site/bctnet/home?authuser=0

load('data/workspace_S132.mat');
% n_cnxs = 306; as per BCT: (each edge is rewired approximately ITER times)
n_iter = 100; % as per Fornito's book chapter 10
n_reals = 1000;

sconnDen_sim = zeros(size(ifod_den_35));

for n = 1:n_reals
    [sconnDen_sim(:,:,n), eff] = randmio_und(ifod_den_35, n_iter);
end

