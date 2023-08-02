classdef count_spk < handle

	methods (Static)

		function [spk_raster,t] = tiny_bin_raster(spikes,edges)
			% count spikes with tiny bins for raster plot
			% Input: spikes - cell array of spike times
			

			% put convolution here as well, set parameters as fixed for now since we're copying chris's settings
			% prepare kernel first
			v = CausalGauss(2.5); % 25 ms for 10 ms bin
			
			% bin spikes with tiny bins for raster
			spk_raster = NaN(numel(spikes),numel(edges)-1);
			for ii = 1:numel(spikes)
				spk_raster(ii,:) = histcounts(spikes{ii},edges); % count spikes
				spk_raster(ii,:) = zscore(spk_raster(ii,:)); % zscore
				spk_raster(ii,:) = conv(spk_raster(ii,:),v,'same'); % convolve
			end
			t = edges(1:end-1) + (edges(2)-edges(1)) / 2;
			

			
			% part of the code that used to go into classifier.plt.posterior_raster
			% [spk_raster,t_raster] = classifier.count_spk.tiny_bin_raster(spk_2_plt,data.video(1):0.01:data.video(end));
			% raster_cmap = flip(cbrewer2('RdBu'));
			
			% % new method, bin 10 ms, convolve with 25 ms half-gaussian - keep constant with Chris method
			% % discarded
			% tmp = (t_raster>toi(1) & t_raster<toi(2));
			% imagesc(ax(ii),t_raster(tmp),1:size(spk_raster,1),spk_raster(:,tmp));
			% colormap(raster_cmap);
			% ax(ii).CLim = [-0.5 0.5];
			% ax(ii).YDir = 'reverse';

		end


		function zscore_by_time = pre_event(data,ops)
			% count spikes for zscoring based on pre-reward or pre-cgrp activities
			
			% partition by time - need the variance from the entire time period
			[~,zscore_by_time] = classifier.count_spk.zscore(data,ops,[common_t.last_reward(data) common_t.first_laser(data)]);

			% calculate baseline for 
			% 1. drinking period - 1 sec precue
			[~,count] = cal_raster(cellfun(@(x) x*1000,data.spikes(~ops.exclude_id),'uni',0), ...
									[data.cues.rewarded.front; data.cues.rewarded.rear]*1000, [-1 0] *1000);
			zscore_by_time.M(1,:) = mean(count,1);
			% 2. 2nd context - keep the original
			% 3. pre-cgrp stim
			[~,count] = cal_raster(cellfun(@(x) x*1000,data.spikes(~ops.exclude_id),'uni',0), ...
									data.laser(:,1)*1000, [-1 0] *1000);
			zscore_by_time.M(3,:) = mean(count,1);


			% plotting the difference in zscoring method
				% [~,zscore_time_course] = classifier.count_spk.zscore(data,ops,[common_t.last_reward(data) common_t.first_laser(data)]);
				% zscore_pre_stim = classifier.count_spk.pre_event(data,ops);
				% plt_cmp_matrices({zscore_time_course.M,zscore_pre_stim.M},{'time course','pre stim'});

		end

		function [spk_count,ops_zscore] = zscore(spk_count,t,t_partition,MV)
			% zscore spk_count matrix by separate time period
			% Input: MV (optional), {M, V} 1*2 cell, use already calculated mean / variance

			% recount data
			% spk_count - data; t - ops
			if isstruct(spk_count) && isstruct(t)
				[spk_count,ops_tmp] = classifier.count_spk.time_course(spk_count,t);
				t = ops_tmp.posterior_t;
			end

			% output zscore mean variance or use input
			if nargin > 3
				ops_zscore.M = MV{1};
				ops_zscore.V = MV{2};
			else
				ops_zscore.M = nan(numel(t_partition)+1,size(spk_count,1));
				ops_zscore.V = ops_zscore.M;
			end
			% wrap partiiton with start and end time
			t_partition = [t(1)-1 t_partition t(end)+1];
			ops_zscore.t_partition = t_partition;


			% do z score
			for ii = 1:(numel(t_partition)-1)
				ind = (t>t_partition(ii))&(t<=t_partition(ii+1));
				if nargin > 3
					spk_count(:,ind) = (spk_count(:,ind)-ops_zscore.M(ii,:)') ./ ops_zscore.V(ii,:)';
				else
					[tmp,ops_zscore.M(ii,:),ops_zscore.V(ii,:)] = zscore(spk_count(:,ind)');
					spk_count(:,ind) = tmp';
				end
			end


		end

		function [count,ops,events_oi] = events(data,ops)
			% count spikes according to events for classifier

			if nargin < 2
				ops = struct;
			end
			ops.tp      = getOr(ops,'tp',[0.1 0.7]);
			ops.events  = getOr(ops,'events',{'front','rear','precue'});; % ops.events  = 
			events_oi   = {data.rewards.all.front,data.rewards.all.rear,[data.cues.rewarded.front; data.cues.rewarded.rear]};

			% make sure only count spikes before cue onset
			t_precue = -flip(ops.tp);
			if t_precue(1,2) > 0
				t_precue = t_precue - (t_precue(1,2));
			end


			% count spikes
			clear count;
			for i = 1:numel(data.spikes)
				for jj = 1:size(ops.tp,1)
					[~,count{1}(jj,i,:)] = cal_raster(data.spikes{i}*1000, events_oi{1}*1000, ops.tp(jj,:) *1000);
					[~,count{2}(jj,i,:)] = cal_raster(data.spikes{i}*1000, events_oi{2}*1000, ops.tp(jj,:) *1000);
					[~,count{3}(jj,i,:)] = cal_raster(data.spikes{i}*1000, events_oi{3}*1000, t_precue *1000);
				end
			end

			% concatenate time windows
			count = arrayfun(@(ii) reshape(count{ii},prod(size(count{ii},[1 2])),[]), 1:3,'UniformOutput',false);

			% remove baseline
			if numel(ops.events)==2
				count(3) = [];
			end
		end


		
		function [spk_count,ops] = time_course(data,ops,recalculate)
			% count spikes for the time course of a session
			% or with edges in ops.posterior_t
			% ops.tp - time window used to count spike for training

			% parameter initiation
			ops.exclude_id    = getOr(ops,'exclude_id',false(size(data.spikes)));
			posterior_t_edges = getOr(ops,'posterior_t_edges',data.video(1):0.5:data.video(end));
			ops.tp = getOr(ops,'tp',[0 1]);
			bin_width = diff(ops.tp); 

			% save temporary file for skipping counting spikes
			tmp = diff(posterior_t_edges(1:2));
			if any(tmp == [0.1 0.5 1]) % need this if because previously only saved 1 decimal
				f_spk = sprintf('mat/bin_%.1f_step_%.1f_%s_%s.mat',diff(ops.tp),diff(posterior_t_edges(1:2)),data.subject,datestr(data.session,'YYmmdd'));
			else
				f_spk = sprintf('mat/bin_%.1f_step_%.2f_%s_%s.mat',diff(ops.tp),diff(posterior_t_edges(1:2)),data.subject,datestr(data.session,'YYmmdd'));
			end

			% load file
			if exist(f_spk) && ((nargin < 3) || ~recalculate)
				load(f_spk);

				% recalculate if not including all neurons
				% if size(spk_count,1)<numel(data.spikes)
				% 	spk_count = NaN(numel(data.spikes),numel(posterior_t_edges)); t = [];
				% 	for ii = 1:numel(data.spikes)
				% 		[~,spk_count(ii,:),~,t] = running_average(data.spikes{ii},[],diff(ops.tp),[],posterior_t_edges); 
				% 	end
				% 	save(f_spk,'spk_count','t','bin_width');
				% end

			% count spikes again
			else
				% use histcount for efficiency
				spk_count = NaN(numel(data.spikes),numel(posterior_t_edges));
				if bin_width == (posterior_t_edges(2)-posterior_t_edges(1))
					for ii = 1:numel(data.spikes)
						[spk_count(ii,:),edges] = histcounts(data.spikes{ii},[posterior_t_edges posterior_t_edges(end)+posterior_t_edges(2)-posterior_t_edges(1)]);
					end
					t = ops.posterior_t_edges + bin_width / 2;
				else
					for ii = 1:numel(data.spikes)
						% [~,spk_count(ii,:),~,t] = running_average(data.spikes{ii},[],diff(ops.tp),[],posterior_t_edges); 
						[spk_count(ii,:),t] = hist_overlap(data.spikes{ii},diff(ops.tp),posterior_t_edges);
					end
				end
				save(f_spk,'spk_count','t','bin_width');
			end

			% select subset
			spk_count = spk_count(~ops.exclude_id,:);
			cell_oi   = find(~ops.exclude_id);

			% zscore by different time periods
			if isfield(ops,'zscore_by_time')
				[spk_count,ops_zscore] = classifier.count_spk.zscore(spk_count,t,ops.zscore_by_time.t_partition(2:end-1),{ops.zscore_by_time.M, ops.zscore_by_time.V});
			end

			% save info
			ops.posterior_t  = t;
			ops.posterior_w  = bin_width; 
			ops.spk_count_id = cell_oi;


		end
	end

end