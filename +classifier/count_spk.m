classdef count_spk < handle

	methods (Static)
		function [count,ops] = events(data,ops)
			% count spikes according to events for classifier

			if nargin < 2
				ops = struct;
			end
			ops.tp      = getOr(ops,'tp',[0.1 0.7]);
			ops.events  = getOr(ops,'events',{'front','rear','precue'});

			% count spikes
			for i = 1:numel(data.spikes)
				[raster,count{1}(i,:)] = cal_raster(data.spikes{i}*1000, data.rewards.all.(ops.events{1})*1000, ops.tp *1000);
				[raster,count{2}(i,:)] = cal_raster(data.spikes{i}*1000, data.rewards.all.(ops.events{2})*1000, ops.tp *1000);
				[raster,count{3}(i,:)] = cal_raster(data.spikes{i}*1000, [data.cues.rewarded.front; data.cues.rewarded.rear]*1000, -flip(ops.tp) *1000);
			end
		end

		function [spk_count,ops] = time_course(data,ops)
			% count spikes for the time course of a session
			% or with edges in ops.posterior_t
			% ops.tp - time window used to count spike for training

			% parameter initiation8
			posterior_t_edges = getOr(ops,'posterior_t_edges',data.video(1):0.5:data.video(end));
			bin_width = diff(ops.tp); 

			spk_count = NaN(sum(~ops.exclude_id),numel(posterior_t_edges)); t = [];
			cell_oi = find(~ops(1).exclude_id);
			for ii = 1:numel(cell_oi)
				[~,spk_count(ii,:),~,t] = running_average(data.spikes{cell_oi(ii)},[],diff(ops.tp),[],posterior_t_edges); 
			end

			% save info
			ops.posterior_t = t;
			ops.posterior_w = bin_width; 


		end
	end

end