

fact_len = getAvgConnectome("Schaefer300_7netANDTianS2_acpc","FACT","SIFT2","length");
fact_den = getAvgConnectome("Schaefer300_7netANDTianS2_acpc","FACT","SIFT2","standard");
fact_len = fact_len([1:166],[1:166]);
fact_den = fact_den([1:166],[1:166]);


% % abagen/nsb have different ordering of aseg
% fact_len_flip = fact_len([1:38,41,39:40],[1:38,41,39:40]);
% fact_den_flip = fact_den([1:38,41,39:40],[1:38,41,39:40]);

fact_den_35 = threshold_proportional(fact_den, 0.35);
fact_len_35 = fact_len;
fact_len_35(fact_den_35 == 0) = 0;
fact_len_35 = reshape(fact_len_35, [166,166]);
% disp(corrcoef(fact_den_35, sconnDen));

fact_den_corr = fact_den_35./fact_len_35;
mask1 = isnan(fact_den_corr);
mask2 = (fact_den_corr == inf);
mask3 = (fact_den_corr == -inf);
fact_den_corr(mask1 | mask2 | mask3) = 0;
% disp(corrcoef(fact_den_corr, sconnDen));

fact_den_corr = (fact_den_corr./ROIsize)'./ROIsize;
mask1 = isnan(fact_den_corr);
mask2 = (fact_den_corr == inf);
mask3 = (fact_den_corr == -inf);
fact_den_corr(mask1 | mask2 | mask3) = 0;
% disp(corrcoef(fact_den_corr, sconnDen));
