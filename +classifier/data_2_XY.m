function [X,Y,ops] = data_2_XY(data,ops)
% function that converts data struct to XY for classifier input
% Output: X - ntrs * ncells spike count
% 		  Y - ntrs * 1 category label

% whether rescale the data based on Andermann et al.
ops.rescale = getOr(ops,'rescale',1);

% convert data to spike count if input is raw data struct
% otherwise it should be the output of count_spk.events
if isstruct(data)
	[spk_count,ops] = classifier.count_spk.events(data,ops);
end

% include certain cells and rescale, default include all cells
% include_id = [9 13 28 29 30 34 46 54 67]; ops.exclude_id = getOr(ops,'exclude_id',~ismember(1:size(count{1},1),include_id));
ops.exclude_id = getOr(ops,'exclude_id',false(1,size(spk_count{1},1)));

% remove cells with 0 variance
% need to move this to naive bayes in the future
if getOr(ops,'if_exclude_0var',false)
	v = catcell(cellfun(@(x) std(x,[],2), spk_count,'UniformOutput',false));
	if size(v,2)~=numel(spk_count) v = v'; end % need this for only 1 cell
	ops.exclude_id = ops.exclude_id | (prod(v,2)==0)';
	ops.exclude_method = [{'0 variance'}, getOr(ops,'exclude_method',{})];
end

% create table and run classifier
data_nb = cellfun(@(x) x(~ops.exclude_id,:) / ops.rescale, spk_count,'UniformOutput',false);
X = catcell(data_nb,2)'; % 90 cells * 100 trials response
Y = catcell(arrayfun(@(i) i*ones(size(data_nb{i},2),1),1:numel(data_nb),'UniformOutput',false),1);


% add fake value to 0 variance dataset
ops.fake0 = getOr(ops,'fake0',0);
for c = unique(Y')
	trs = find(Y==c); % find all samples within a category
	zero_ind = find(std(X(trs,:),[],1)==0); % cells with 0 variance
	X(trs(end),zero_ind) = X(trs(end),zero_ind) + ops.fake0; % add 1e-5 to the last trial within the category
end

% keep final neurons in the decoder
ops.decoder_id = find(~ops.exclude_id);

end