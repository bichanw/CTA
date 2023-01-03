function [count,ops] = count_spk(data,ops)
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

