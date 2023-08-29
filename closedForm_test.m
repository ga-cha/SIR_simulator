% main.m
% a script to simulate atrophy accrual due to the accumulation of misfolded alpha-syn aggregates
% 42-region parcellation
% load gene expressions, real atrophy, ROIsize, functional connectivity...
load('data_yqz/42regions/workspace.mat');
% load structural connectivity
load('data_yqz/42regions/sc30.mat');

N_regions = 42;
v = 1;
dt = 0.01;
T_total = 20000;
init_number = 1;
syn_control = ROIsize;
prob_stay = 0.5;
trans_rate = 1;
seed = N_regions;
% load your GBA, SNCA, sconnDen, sconnLen, ROISize ....

allclose = @(x, y, tol) all(abs(x-y)<tol, 'all');
imageFunction = @(a,b) eval('figure; subplot(1, 3, 1); imagesc(a); colorbar; subplot(1,3,2); imagesc(b); colorbar; subplot(1, 3, 3); imagesc(a-b); colorbar');

%%
%%%%% simulation ------ >>>

f1 = @() SIRsimulator(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate);
f2 = @() SIRsimulator3(N_regions, v, dt, T_total, GBA, SNCA, sconnLen, sconnDen, ROIsize, seed, syn_control, init_number, prob_stay, trans_rate);

clc
tic
[Rnor_all, Rmis_all, Rnor0] = f1();
toc

tic
[Rnor_all_2, Rmis_all_2, Rnor0_2] = f2();
toc

% timeit(f1), timeit(f2)

allclose(Rnor0, Rnor0_2, 1e-7)
allclose(Rnor_all, Rnor_all_2, 1e-7)
allclose(Rmis_all, Rmis_all_2, 1e-7)


imageFunction(Rnor_all, Rnor_all_2);
imageFunction(Rmis_all, Rmis_all_2);



%%
figure; 
subplot(2, 4, 1); imagesc(Rnor_all'); title('R_{nor}');
subplot(2, 4, 2); imagesc(Rnor_all_2'); title('R_{nor}^{test}');
subplot(2, 4, 3); imagesc(Rmis_all'); title('R_{mis}')
subplot(2, 4, 4); imagesc(Rmis_all_2'); title('R_{mis}^{test}');

subplot(2, 4, [5 6]); imagesc(abs(Rnor_all' - Rnor_all_2')./Rnor_all'<0.1); 
title("Where is R_{nor}^{test} within 10% of expected?")

subplot(2, 4, [7 8]); imagesc(abs(Rmis_all' - Rmis_all_2')./Rmis_all'<0.1); 
title("Where does R_{mis}^{test} within 10% of expected?")

%%
ratio = Rmis_all ./(Rnor_all + Rmis_all) ;
ratio(Rmis_all<1) = 0; % remove possible NaNs...

% atrophy growth
k1 = 0.5;
k2 = 1 - k1;
% input weigths of deafferentation (scaled by structrual connectivity)
weights = sconnDen ./ repmat(sum(sconnDen, 2), 1, N_regions);

% neuronal loss caused by lack of input from neighbouring regions
ratio_cum = weights * (1-exp(-ratio * dt));
% one time step back
ratio_cum = [zeros(N_regions, 1), ratio_cum(:, 1:end-1)];
ratio_cum = k2 * ratio_cum + k1 * (1-exp(-ratio * dt));

% add all the increments across t
simulated_atrophy = cumsum(ratio_cum, 2);
