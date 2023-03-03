classdef count_spk < handle

	methods (Static)
		function [count,ops,events_oi] = events(data,ops)
			% count spikes according to events for classifier

			if nargin < 2
				ops = struct;
			end
			ops.tp      = getOr(ops,'tp',[0.1 0.7]);
			ops.events  = {'front','rear','precue'}; % ops.events  = getOr(ops,'events',{'front','rear','precue'});
			events_oi   = {data.rewards.all.front,data.rewards.all.rear,[data.cues.rewarded.front; data.cues.rewarded.rear]};

			% count spikes
			clear count;
			for i = 1:numel(data.spikes)
				for jj = 1:size(ops.tp,1)
					[~,count{1}(jj,i,:)] = cal_raster(data.spikes{i}*1000, events_oi{1}*1000, ops.tp(jj,:) *1000);
					[~,count{2}(jj,i,:)] = cal_raster(data.spikes{i}*1000, events_oi{2}*1000, ops.tp(jj,:) *1000);
					[~,count{3}(jj,i,:)] = cal_raster(data.spikes{i}*1000, events_oi{3}*1000, -flip(ops.tp(jj,:)) *1000);
				end
			end

			% concatenate time windows
			count = arrayfun(@(ii) reshape(count{ii},prod(size(count{ii},[1 2])),[]), 1:3,'UniformOutput',false);
			% count{1}(i,j,k) -> tmp(2*(j-1)+i)
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
			f_spk = sprintf('mat/bin_%.1f_step_%.1f_%s_%s.mat',diff(ops.tp),diff(posterior_t_edges(1:2)),data.subject,datestr(data.session,'YYmmdd'));

			% load file
			if exist(f_spk) || (nargin>2 && ~recalculate)
				load(f_spk);
			% count spike
			else
				spk_count = NaN(sum(~ops.exclude_id),numel(posterior_t_edges)); t = [];
				cell_oi = find(~ops.exclude_id);
				for ii = 1:numel(cell_oi)
					[~,spk_count(ii,:),~,t] = running_average(data.spikes{cell_oi(ii)},[],diff(ops.tp),[],posterior_t_edges); 
				end

				save(f_spk,'spk_count','t','bin_width','cell_oi');
			end

			% save info
			ops.posterior_t = t;
			ops.posterior_w  = bin_width; 
			ops.spk_count_id = cell_oi;


		end
	end

end