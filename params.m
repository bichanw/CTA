ops    = struct('if_cv',false,'classifier',classifier.mnr(),'tp',[0 1]);  % classifier
prefix = '';

% multinomial logistic regression settings
ops.mnr.lambda = 1;
ops.mnr.penalty = 'l1';

% select significant cells
ops.novel_vs_fam.n_sig = 15;
ops = classifier.select_cells.novel_vs_fam(data,ops);
% ops.exclude_id = ~ismember(1:numel(data.spikes),ops.novel_vs_fam.ordered_id(1:ops.novel_vs_fam.ordered_div(end-2)));
% ops.exclude_method = {'novel_vs_fam'};
ops = classifier.select_cells.sig_resp(data,ops); % select all significant


% some plotting options
ops.plot = struct();
ops.plot.slow_firing_cell = 'novel_vs_fam'; %{'novel_vs_fam','non_zero_coef'}


% do zscore based on separate time scale
% [~,ops.zscore_by_time] = classifier.count_spk.zscore(data,ops,[common_t.last_reward(data) common_t.first_laser(data)]);
% ops = classifier.select_cells.non_zero_zscore(ops); % clean up cells with 0 variance
ops.mnr.zscore = true; % do not do separate zscore


% remove baseline in the classifier
% ops.events = {'front','rear'};


% set of available options
% classifier: object of a classifier with defined train / predict methods
% exclude_id: boolean vector of n cells (same in data.spikes) * 1, which cells to exclude
% 
% novel_vs_fam: struct containing options / results for the ranksum test
% 				n_sig: 	 	choose the n most signifcant cells
% 				ordered_id: the index of cells (i in data.spikes{i}) that pass the ranksum test 
% 				
% decoder_id: defined in data_2_XY.m, data.spikes(ops.decoder_id) was used in the decoder analysis
% 
% by_brain_region: {'BLA','CEA'}
% 
% z_score_by_time: mean / variance for zscoring across time partition
% 				   fields:   M, V, t_partition

