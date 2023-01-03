function [cell_i,data] = sort_raster(pref,data,pref_ind,count)

% find significant cells
cids = pref.clusterid;
cell_i = arrayfun(@(x) find(data.cids==x), cids);

% sort cells for raster
sorting_method = 'dp';
switch sorting_method
	case 'acc'
		[~,I] = sort(pref.accuracy);
		cell_i = cell_i(I);
	case 'dp'
		dp = my_dp(count{pref_ind}',count{3}');
		% [~,I] = sort(abs(dp(cell_i)));
		[~,I] = sort(dp(cell_i));
		cell_i = cell_i(I(1:min([10 numel(I)]))); % keep at most 10 cells
		% cell_i = cell_i(I); % keep all cells
end

% keep group info
data.group = [data.group; pref_ind * ones(numel(cell_i),1)];


% dp = my_dp(count{1}',count{3}');