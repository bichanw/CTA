classdef select_cells < handle
% master class for naive bayes classifier


	% by cell
	methods (Static)

		function ops = novel_vs_fam(data,ops)
			% count spikes 
			[spk_count,ops] = classifier.count_spk.events(data,ops);

			% categorize cell
			n_cells  = size(spk_count{1},1);
			cell_cat = nan(n_cells,1); % cell category - 0 - insignificant
			zval     = cell_cat;
			for ii = 1:n_cells
				% novel vs familiar
				[p,is_sig,stats] = ranksum(spk_count{1}(ii,:),spk_count{2}(ii,:));


				% if response is different
				if is_sig
					% front activation higher: stats.zval > 0, category-1
					% back activation higher: stats.zval < 0, category-2
					cell_cat(ii) = (stats.zval<0) + 1;
					zval(ii)     = stats.zval;

				% no difference between front vs back
				else
					[p,is_sig,stats] = ranksum([spk_count{1}(ii,:) spk_count{2}(ii,:)],spk_count{3}(ii,:));
					if is_sig 
						cell_cat(ii) = 3; % different than control
						zval(ii)     = stats.zval;
					else
						cell_cat(ii) = 0; % no response
					end

				end
			end

			% pick n most signifcant ones
			ops.novel_vs_fam = getOr(ops,'novel_vs_fam',struct());
			n_sig = getOr(ops.novel_vs_fam,'n_sig',10);
			ordered_id  = [];
			ordered_div = 0; % category divider
			for ii = 1:3 % 3 category of interest
				% pick out cells
				ind = find(cell_cat==ii);

				% rank cell by absolute zval
				[B,I] = sort(abs(zval(ind)),'descend');
				n = min([numel(I) n_sig]);
				ordered_id = [ordered_id; ind(I(1:n))];

				% keep divider info
				ordered_div = [ordered_div ordered_div(end)+n];
			end

			% add non significant cells
			ind = find(cell_cat==0);
			n = min([numel(ind) n_sig]);
			ordered_id = [ordered_id; ind(randperm(numel(ind),n))];
			ordered_div = [ordered_div ordered_div(end)+n];

			% save information
			ops.novel_vs_fam.ordered_id = ordered_id;
			ops.novel_vs_fam.ordered_div = ordered_div;
			ops.novel_vs_fam.n_sig = n_sig;
			ops.novel_vs_fam.cell_cat = cell_cat;
			ops.novel_vs_fam.cell_cat_name = {'front higher','back higher','different from control','non significant'};
			
			
		end

		function ops = sig_resp(data,ops)
			% select cells with significant response to reward

			% count spikes 
			[spk_count,ops] = classifier.count_spk.events(data,ops);

			% test for significance
			n_cells = size(spk_count{1},1);
			pairs   = [1 3;2 3;1 2];
			is_sig  = nan(n_cells,size(pairs,1));
			for ii = 1:n_cells
				for jj = 1:size(pairs,1)
					[~,is_sig(ii,jj)] = ranksum(spk_count{pairs(jj,1)}(ii,:),spk_count{pairs(jj,2)}(ii,:));
				end
			end


			% exclude cells with none significant
			ops.exclude_id = getOr(ops,'exclude_id',false(1,numel(data.spikes))) | ~any(is_sig,2)';
			ops.exclude_method = [{'significant response'}, getOr(ops,'exclude_method',{})];


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