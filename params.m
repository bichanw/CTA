ops    = struct('if_cv',false,'classifier',classifier.mnr(),'tp',[0 1]);  % classifier
prefix = '';

% multinomial logistic regression settings
ops.mnr.lambda = 1;
ops.mnr.penalty = 'l1';
ops.amplitude_cutoff = 20;
ops.amplitude_cat = {'novel','water','neither'}; % categories of clusters to include: {'novel','water','neither'}


% amplitude cutoff -- need to be at the end of script after other criterion screening!!!
ops.amplitude_cutoff = getOr(ops,'amplitude_cutoff',20);
if ops.amplitude_cutoff >= 0 % amplitude cutoff
	ops.exclude_method = {'amplitude'};
	load(sprintf('/jukebox/witten/Chris/matlab/cz/neuropixels-cta/bichan/calca%s_clusters_%duV_FDR5pct_10sec.mat',data.subject,ops.amplitude_cutoff))
	include_id = [];
	for icat = 1:numel(ops.amplitude_cat)
		include_id = [include_id, cluster_assignments.(ops.amplitude_cat{icat})];
	end
	ops.exclude_id = ~ismember(data.cids,include_id)';
	% ops = classifier.select_cells.sig_resp(data,ops); % select all significant - did not change results by much

else % or select cells with significant response
	ops = classifier.select_cells.sig_resp(data,ops); % select all significant
end

% 
if strcmpi(data.subject,'302')
	ops.exclude_id(100) = true;
end


% or selecting cells that is most significant
% ops.exclude_id = ~ismember(1:numel(data.spikes),ops.novel_vs_fam.ordered_id(1:ops.novel_vs_fam.ordered_div(end-2)));

% categorizing cells for plotting purposes
ops.novel_vs_fam.n_sig = 15;
ops = classifier.select_cells.novel_vs_fam(data,ops);
	
% some plotting options
ops.plot = struct();
ops.plot.slow_firing_cell = 'novel_vs_fam'; %{'novel_vs_fam','non_zero_coef'}


% zscore
ops.zscore_method = '2nd_context';
switch ops.zscore_method
case 'separate' % do zscore based on separate time scale
	% [spk_count_zscored,ops.zscore_by_time] = classifier.count_spk.zscore(data,ops,[common_t.last_reward(data) common_t.first_laser(data)]);
	[~,ops.zscore_by_time] = classifier.count_spk.zscore(data,ops,[common_t.last_reward(data) common_t.first_laser(data)]);
	ops = classifier.select_cells.non_zero_zscore(ops); % clean up cells with 0 variance
	ops.mnr.zscore = false;
case 'all' % zscore based on the same mean / variance
	ops.mnr.zscore = true; 
case 'pre_event'
	ops.zscore_by_time = classifier.count_spk.pre_event(data,ops);
	ops.mnr.zscore = false;
case '2nd_context'
	% use separate period, but using 2nd context for both 2nd context and stim period
	ops.zscore_by_time = classifier.count_spk.pre_event(data,ops);
	ops.zscore_by_time.M(3,:) = ops.zscore_by_time.M(2,:);
	ops.zscore_by_time.V(3,:) = ops.zscore_by_time.V(2,:);
	ops.mnr.zscore = false;
end





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

