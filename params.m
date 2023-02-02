ops = struct('if_cv',false,'classifier',classifier.mnr());  % classifier

% multinomial logistic regression settings
ops.mnr.lambda = 1;
ops.mnr.penalty = 'l2';


% set of available options
% classifier: object of a classifier with defined train / predict methods
% exclude_id: boolean vector of n cells (same in data.spikes) * 1, which cells to exclude
% 
% novel_vs_fam: struct containing options / results for the ranksum test
% 				n_sig: 	 	choose the n most signifcant cells
% 				ordered_id: the index of cells (i in data.spikes{i}) that pass the ranksum test 
% 				