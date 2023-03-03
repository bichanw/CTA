ops = struct('if_cv',false,'classifier',classifier.mnr(),'tp',[0 1]);  % classifier

% multinomial logistic regression settings
ops.mnr.lambda = 1;
ops.mnr.penalty = 'l1';


% select significant cells
ops.novel_vs_fam.n_sig = 15;
ops = classifier.select_cells.novel_vs_fam(data,ops);
% ops.exclude_id = ~ismember(1:numel(data.spikes),ops.novel_vs_fam.ordered_id(1:ops.novel_vs_fam.ordered_div(end-2)));
% ops.exclude_method = {'novel_vs_fam'};
ops = classifier.select_cells.sig_resp(data,ops); % select all significant


% set of available options
% classifier: object of a classifier with defined train / predict methods
% exclude_id: boolean vector of n cells (same in data.spikes) * 1, which cells to exclude
% 
% novel_vs_fam: struct containing options / results for the ranksum test
% 				n_sig: 	 	choose the n most signifcant cells
% 				ordered_id: the index of cells (i in data.spikes{i}) that pass the ranksum test 
% 				
% decoder_id: defined in data_2_XY.m, data.spikes(ops.decoder_id) was used in the decoder analysis