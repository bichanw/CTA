classdef select_cells < handle
% master class for naive bayes classifier


	% by cell
	methods (Static)

		function ops = sig_resp(data,ops)
			% select cells with significant response to reward

			% count spikes 
			[spk_count,ops] = classifier.count_spk.events(data,ops);

			% test for significance
			n_cells = size(spk_count{1},1);
			is_sig  = nan(n_cells,2);
			for ii = 1:n_cells
				for jj = 1:2
					[p,is_sig(ii,jj)] = ranksum(spk_count{jj}(ii,:),spk_count{3}(ii,:));
				end
			end

			% exclude cells with neither significant
			ops.exclude_id = getOr(ops,'exclude_id',false(1,numel(data.spikes)));
			ops.exclude_id = ops.exclude_id | ~(is_sig(:,1)|is_sig(:,2))';

		end

		function ops = individual_cv(data,ops)

			% parameter
			thresh = 0.6;
			ops.select_cell = sprintf('misclassification < %.2f',thresh);

			% run classifier with individual cell and exclude bad ones
			[cost,valid] = classifier.nb.individual_cv(data,ops);
			ops.exclude_id = getOr(ops,'exclude_id',false(1,numel(data.spikes)));
			ops.exclude_id = ops.exclude_id | ~(cost'<thresh);
		end
	end
	

end