function [R,count,edges,tbin,V] = running_average(T,data,bin_width,step_size,edges,otherfunc)
% [R,count,edges,tbin] = running_average(T,data,bin_width,step_size,edges)
% calculate the running average of data
% each data happens at T
% each bin has count(ibin) number of data points
if isempty(data)
	data = ones(size(T));
end
if nargin < 5 || isempty(edges)
	edges = [min(T):step_size:(max(T)-bin_width)];
elseif isempty(bin_width)
	bin_width = edges(2) - edges(1);
end

NBins = numel(edges);

% let the first dim of T be NBins
if size(data,1) ~= numel(T)
	data = data';
end
NDims = size(data,2);

R = nan(NBins,NDims);
count = zeros(NBins,1);

for iBin = 1:NBins
	% find data in the bin
	data_in_bin = T>=edges(iBin)&T<(edges(iBin)+bin_width);

	% print message if data points too few
	count(iBin) = sum(data_in_bin);
	if sum(data_in_bin)<5
		% fprintf('Few than 5 data in bin %d\n',iBin);

		% skip no data point
		% if sum(data_in_bin) < 1
			continue;
		% end
	end

	% calculate average
	if nargin < 6
		R(iBin,:) = nanmean(data(data_in_bin,:),1);
		V(iBin,:) = std(data(data_in_bin,:),[],1) / count(iBin);
	else
		R(iBin,:) = otherfunc(data(data_in_bin,:));
	end
end

if nargout > 3
	tbin = edges + bin_width /2;
end