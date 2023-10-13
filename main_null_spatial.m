% main_null_spatial.m
%                                                                                                                                    
% Package based off https://github.com/yingqiuz/SIR_simulator
% Zheng, Ying-Qiu, et al. PLoS biol. 17.11 (2019): e3000495.
%
% Gabriella Chan 12/10/23
% gabriella.chan@monash.edu
% Monash University
%
% We generate the simulated atrophy map for a gene pair, and compare with
% spatial nulls. Nulls are scrambles of empirical atrophy, simulated
% with BrainSMASH (brainsmash.readthedocs.io)
% We calculate the p-value of empirical atrophy against the null, and plot
% on a box-and-whiskers plot, with whiskers set to the 95% CI

% load gene expressions, real atrophy, ROIsize, functional connectivity...
load('data_gc/GC_workspace.mat');

N_regions = 41;
v = 1;
dt = 0.01;
T_total = 10000;
init_number = 1;
syn_control = ROIsize;
prob_stay = 0.5;
trans_rate = 1;
% init seed to hip
seed = 40;
% single risk/clearance gene pair input, as tables
clear_gene = genes(:, 'REEP4');
risk_gene = LAMP5;

% First we generate simulated atrophy 
[Rnor_all, Rmis_all] = SIRsimulator(N_regions, v, dt, T_total, ...
    clear_gene, risk_gene, sconnLen, sconnDen, ROIsize, seed, ...
    syn_control, init_number, prob_stay, trans_rate);
[sim_atrophy] = SIRatrophy(Rnor_all, Rmis_all, sconnDen, N_regions, dt);

% Then we determine the correlation with null atrophy. Null atrophy is an 
% input of scrambled empirical atrophy, preserving spatial autocorrelation

nulls = 1000;
null_corrs = zeros(nulls,1);

for i = 1:nulls
    [null_corrs(i), ~] = SIRcorr(sim_atrophy, null_atrophy(i, :), T_total);
end

% experimental atrophy correlation
[exp_corr, ~] = SIRcorr(sim_atrophy, emp_atrophy, T_total);

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
boxplot(null_corrs, 'whisker', w95)
swarmchart(ones(length(null_corrs),1),null_corrs,10,'jitter','on')
ylim([0 0.85])
scatter(1,exp_corr,20,"filled")
hold off