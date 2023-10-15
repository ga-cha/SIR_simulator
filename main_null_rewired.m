% main_null_rewired.m
%                                                                                                                                    
% Package based off https://github.com/yingqiuz/SIR_simulator
% Zheng, Ying-Qiu, et al. PLoS biol. 17.11 (2019): e3000495.
%
% Gabriella Chan 13/08/23
% gabriella.chan@monash.edu
% Monash University
%
% We generate simulated atrophy of our gene pair with each different
% rewired null, and record the correlation with empirical atrophy. We then
% compute the correlation of simulated atrophy with empirical atrophy when
% using the true structural connectome density.
% We calculate the p-value of experimental correlation against nulls, and 
% plot on a box-and-whiskers plot.

% load gene expressions, real atrophy, ROIsize, functional connectivity...
load('data_gc/GC_workspace.mat');

% Load simulated structural connectivity, since Maslov-Sneppen may connect
% unconnected nodes. Described in greater detail in Zheng PLoS biol. (2019)
% The implementation is given in SIR_utils/null_ROI_dist.m
load('data_gc/sconnLen_sim.mat');

N_regions = 41;
v = 1;
dt = 0.01;
T_total = 10000;
init_number = 1;
syn_control = ROIsize;
prob_stay = 0.5;
trans_rate = 1;
% initialise seed to hip
seed = 40;
% single risk/clearance gene pair input, as tables
clear_gene = genes(:, 'GNLY');
risk_gene = genes(:, 'EHMT2');
nulls = 1000;
null_corrs = zeros(nulls,1);

% initialise the array
for i = 1:nulls
    f = strcat('data_gc/rewire/rewire', num2str(i), '.csv');
    sconnDen_sim = readmatrix(f);
    sconnDen_sim = sconnDen_sim(1:N_regions,1:N_regions);
    [gene_corrs, ~] = SIRiterator(N_regions, v, dt, T_total, ...
        clear_gene, risk_gene, sconnLen_sim, sconnDen_sim, ROIsize, seed, ...
        syn_control, init_number, prob_stay, trans_rate, emp_atrophy);
    null_corrs(i) = gene_corrs.correlation;
end

% experimental atrophy
[gene_corrs, ~] = SIRiterator(N_regions, v, dt, T_total, ...
    clear_gene, risk_gene, sconnLen, sconnDen, ROIsize, seed, ...
    syn_control, init_number, prob_stay, trans_rate, emp_atrophy);
exp_corr = gene_corrs.correlation;
% p-value
[h,p]=ttest2(exp_corr,null_corrs);
disp(['p = ', num2str(p)])

% Set box plot whiskers to 95% CI. Taken from:
% https://www.mathworks.com/matlabcentral/answers/171414-how-to-show-95-quanile-in-a-boxplot
q95=norminv(0.95);
q3=norminv(.75);
w95=(q95-q3)/(2*q3);

figure;
hold on
boxplot(null_corrs, 'whisker', 0.7193)
swarmchart(ones(length(null_corrs),1),null_corrs,10,'jitter','on')
ylim([0 0.85])
scatter(1,exp_corr,20,"filled")
hold off