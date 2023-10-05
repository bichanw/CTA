classdef select_cells < handle
% master class for naive bayes classifier


	% by cell
	methods (Static)


		function ops = non_zero_zscore(ops)
			ind_all = find(~ops.exclude_id);
			ind_exlcude = find(logical(sum(ops.zscore_by_time.V==0,1)));

			% exclude cell
			ops.exclude_id(ind_all(ind_exlcude)) = true;
			% clean up 
			ops.zscore_by_time.M(:,ind_exlcude) = [];
			ops.zscore_by_time.V(:,ind_exlcude) = [];
		end

		function ind = by_brain_region(data,ops)
			% Input:  ops - either ops struct of cells of regions of interest
			% Output: ind - either cell (multiple roi) or vector (1 roi)

			% extract unit's corresponding region
			units_region = data.brain_region(sub2ind(size(data.brain_region),data.sh+1,data.ch+1));
		
			% find regions
			if isstruct(ops)
				regions_oi = getOr(ops,'by_brain_region',{'CEA'});
			else
				regions_oi = ops;
			end

			% all units in region of interest
			for ii = 1:numel(regions_oi)
				ind{ii} = cellfun(@(x) contains(x,regions_oi{ii}), units_region);
			end

			% compress for easier processing
			if numel(ind)==1  ind = ind{1};  end

		end

		function cell_order = by_chris(data,ops)
			% categorize cells by chris's method
			load(sprintf('/jukebox/witten/Chris/matlab/cz/neuropixels-cta/bichan/calca%s_clusters_%duV_FDR5pct_10sec.mat',data.subject,ops.amplitude_cutoff));
			% fprintf('/jukebox/witten/Chris/matlab/cz/neuropixels-cta/bichan/calca%s_clusters_%duV_FDR5pct_10sec.mat',data.subject,ops.amplitude_cutoff);
			% find included cells
			ind  = find(~ops.exclude_id);

			% find order within each category
			ordered_id  = [];
			ordered_div = 0;
			% ordered_id needs to be front vs back
			events_oi   = {data.rewards.all.front,data.rewards.all.rear,[data.cues.rewarded.front; data.cues.rewarded.rear]};
			cells_oi = arrayfun(@(ii) cluster_assignments.(ops.amplitude_cat{ii}), 1:numel(ops.amplitude_cat),'uni',0);
			cells_oi(1:2) = cells_oi(data.port_is_water+1);

			for ii = 1:numel(ops.amplitude_cat)
				% find cells
				tmp = find(ismember(data.cids,cluster_assignments.(ops.amplitude_cat{ii})));
				% rank by peak skipped for the moment?
				if ii < 3
					tmp = tmp(classifier.select_cells.rank_by_peak(data.spikes(tmp),events_oi{ii},ops));
				else
					% subselect control cells
					FR_response = cal_FR(data.spikes(ordered_id));
					ind_2_plt = tmp(randperm(numel(tmp),15)); % select 10 control cells
					FR_control = cal_FR(data.spikes(ind_2_plt));
					
					% make sure the firing rates of control cells match the responsive ones
					while ttest2(FR_response,FR_control)
						ind_2_plt = tmp(randperm(numel(tmp),15)); % select 10 control cells
						FR_control = cal_FR(data.spikes(ind_2_plt));
					end
					tmp = ind_2_plt;
					% ax = np(2,1); histogram(ax(1),FR_response,0:1:30); histogram(ax(2),FR_control,0:1:30); title(ax(1),'responsive');title(ax(2),'control');export_fig('ttest2.pdf');
				end
				ordered_id  = [ordered_id; tmp];
				ordered_div = [ordered_div ordered_div(end)+numel(tmp)];
			end

			% output
			cell_order.ordered_id  = ordered_id;
			cell_order.ordered_div = ordered_div;
			cell_order.cell_cat_name = {'front','back','control'};
		end

		function cell_order = by_rastermap(data,ops)
			% sort cells unsupervised by rastermap

			% save spike count during drinking period for rastermap
			t = (common_t.first_reward(data)-5):0.1:(common_t.last_reward(data)+5);
			spks = nan(sum(~ops.exclude_id),numel(t)-1);
			ind  = find(~ops.exclude_id);
			for ii = 1:size(spks,1)
				spks(ii,:) = histcounts(data.spikes{ind(ii)},t);
			end
			writeNPY(spks,'mat/spks.npy');

			% run rastermap sort
			system('python helpfun/sort_rastermap.py');
			% read rastermap result
			cell_order.ordered_id = ind(readNPY('mat/sorted_rastermap.npy')+1);
			cell_order.ordered_div = [0 numel(cell_order.ordered_id)];
			cell_order.cell_cat_name = {'all'};

		end

		function cell_order = non_zero_coef(data,ops)
			% run classifier
			if ~isfield(ops,'Mdl')
				[~,ops] = ops.classifier.train(data,ops);
			end
			ndim = size(ops.Mdl.coef,2);
			if ndim>1
				d_coef = ops.Mdl.coef(:,1) - ops.Mdl.coef(:,2);
			else
				d_coef = -ops.Mdl.coef;
			end

			cell_order.ordered_id  = [];
			cell_order.ordered_div = 0;
			for ii = 1:2
				% find cells
				if ndim > 1
					ind = find(ops.Mdl.coef(:,ii)~=0);
				else
					ind = find(d_coef * (-2*ii+3) > 0);
				end

				% rank by d_coef (reverse for rear preferred)
				[~,I] = sort(d_coef(ind) * (-2*ii+3),'descend');

				cell_order.ordered_id  = [cell_order.ordered_id ops.decoder_id(ind(I))];
				cell_order.ordered_div = [cell_order.ordered_div numel(cell_order.ordered_id)];
			end

			% replace results
			cell_order.cell_cat_name = {'front nonzero','back nonzero'};

		end

		function [ops,cell_order] = by_coef(data,ops)
			% sort cells by classifier coeffieicent

			% run classifier
			if ~isfield(ops,'Mdl')
				[~,ops] = ops.classifier.train(data,ops);
			end

			% front vs back
			d_coef = ops.Mdl.coef(:,1) - ops.Mdl.coef(:,2);
			[~,I] = sort(d_coef,'descend');
			ind   = find(~ops.exclude_id);

			% pick out cells that preferred novel or familiar based on classifier
			n_sig = 10;
			cell_oi = {ind(I(1:n_sig)),ind(I(end-(n_sig-1):end))};
			cell_oi{1} = cell_oi{1}(ops.Mdl.coef(I(1:n_sig),1)>0);
			cell_oi{2} = cell_oi{2}(ops.Mdl.coef(I(end-(n_sig-1):end),2)>0);

			% sort by peak
			[~,ops,events_oi] = classifier.count_spk.events(data,ops);
			ind_sort_by_peak = [];
			for ii = 1:2
				tmp = cell_oi{ii};
				tmp = tmp(classifier.select_cells.rank_by_peak(data.spikes(tmp),events_oi{ii},ops));
				ind_sort_by_peak = [ind_sort_by_peak, tmp];
			end
			% add group of random cells
			[~,I] = sort(abs(d_coef),'ascend');
			ind_sort_by_peak = [ind_sort_by_peak, ind(I(1:n_sig))];
			

			% replace results
			ops.novel_vs_fam.n_sig         = n_sig;
			ops.novel_vs_fam.ordered_id    = ind_sort_by_peak;
			ops.novel_vs_fam.ordered_div   = [0 cumsum([numel(cell_oi{1}) numel(cell_oi{2}) n_sig])];
			ops.novel_vs_fam.cell_cat_name = {'front higher','back higher','control'};
			if nargout > 1
				cell_order = ops.novel_vs_fam;
			end

		end

		function [ops,cell_order] = novel_vs_fam(data,ops)
			% count spikes 
			[spk_count,ops,events_oi] = classifier.count_spk.events(data,ops);
			ops.exclude_id = getOr(ops,'exclude_id',false(size(spk_count{1},1),1));

			% categorize cell based on rank sum
			n_cells  = size(spk_count{1},1);
			cell_cat = nan(n_cells,1); % cell category - 0 - insignificant
			zval     = cell_cat;
			for ii = 1:n_cells

				% check if already excluded
				if ops.exclude_id(ii)
					cell_cat(ii) = -1;
					continue;
				end

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
						cell_cat(ii) = 3; % no response
					end

				end
			end

			% pick n most signifcant ones
			ops.novel_vs_fam = getOr(ops,'novel_vs_fam',struct());
			n_sig = getOr(ops.novel_vs_fam,'n_sig',15);
			ordered_id  = [];
			ordered_div = 0; % category divider
			include_id  = find(~getOr(ops,'exclude_id',false(numel(data.spikes),1)));
			for ii = 1:3 % 3 category of interest
				% pick out cells
				ind = find(cell_cat==ii);


				% pick most significant cells
				[B,I] = sort(abs(zval(ind)),'descend');
				n = min([numel(I) n_sig]);
				tmp = ind(I(1:n));

				% then rank by peak
				if ii < 3
					tmp = tmp(classifier.select_cells.rank_by_peak(data.spikes(tmp),events_oi{ii},ops));
				end
				tmp = tmp(ismember(tmp,include_id));

				ordered_id = [ordered_id; tmp];
				% keep divider info
				ordered_div = [ordered_div ordered_div(end)+numel(tmp)];
			end

			% add non significant cells
			% removed, since we want to combine non-responsive cells with non-different cells
			% ind = find(cell_cat==0);
			% n = min([numel(ind) n_sig]);
			% ordered_id = [ordered_id; ind(randperm(numel(ind),n))];
			% ordered_div = [ordered_div ordered_div(end)+n];

			% save information
			ops.novel_vs_fam.ordered_id = ordered_id;
			ops.novel_vs_fam.ordered_div = ordered_div;
			ops.novel_vs_fam.n_sig = n_sig;
			ops.novel_vs_fam.cell_cat = cell_cat;
			ops.novel_vs_fam.cell_cat_name = {'front higher','back higher','control'};
			
			if nargout > 1
				cell_order = ops.novel_vs_fam;
			end
			
		end


		function I = rank_by_peak(Spiketimes,eventtime,ops)
			% Spiketimes - nneuron * 1 cell, spike time of each neuron
			% eventtime  - vector of event times
			% ops        - settings inherited 

			% initiation
			tp = getOr(ops,'tp',[0 1]);
			tp = tp(1):0.01:tp(end);

			% calculate psth and order
			resp = catcell(cellfun(@(x) cal_psth(x*1000,eventtime*1000,'tp',tp,'kernel_width',0.1),Spiketimes,'UniformOutput',false));
			[~,max_ti] = max(abs(resp),[],2);
			[~,I] = sort(max_ti,'ascend');
		end

		function ops = sig_resp(data,ops)
			% select cells with significant response to reward

			% count spikes 
			[spk_count,ops] = classifier.count_spk.events(data,ops);

			% test for significance
			n_cells = size(spk_count{1},1);
			% pairs   = [1 3;2 3;1 2];
			pairs   = [1 2];
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