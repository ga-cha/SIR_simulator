% null_rewire
%
% Gabriella Chan 01/08/23
% gabriella.chan@monash.edu
% Monash University
%
% This module calls randmio_und.m from the Brain Connectivity Toolbox to
% generate rewired null networks using the Maslov-Sneppen algorithm
% https://sites.google.com/site/bctnet/home?authuser=0

load('data_gc/GC_workspace.mat');
n_cnxs = 306;
n_iter = n_cnxs * 100; % as per Fornito's book chapter 10
n_reals = 10000;

for n = 4850:n_reals
    [R, S, eff] = randmio_und_den_len(sconnDen, sconnLen, n_iter);
    writematrix(R, strcat('data_gc/rewire/rewireDen', num2str(n), '.csv'));
    writematrix(S, strcat('data_gc/rewire/rewireLen', num2str(n), '.csv'));
end

